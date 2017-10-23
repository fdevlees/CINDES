use_gpu=False
import numpy as np
from sklearn.kernel_ridge import KernelRidge

from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.utils import check_X_y
from sklearn.utils.validation import check_is_fitted

from scipy import linalg

#from sklearn.linear_model.ridge import _solve_cholesky_kernel

class GPU_KernelRidge(KernelRidge):

    def fit(self, X, y=None, sample_weight=None):
        """Fit Kernel Ridge regression model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : array-like, shape = [n_samples] or [n_samples, n_targets]
            Target values

        sample_weight : float or numpy array of shape [n_samples]
            Individual weights for each sample, ignored if None is passed.

        Returns
        -------
        self : returns an instance of self.
        """
        # Convert data
        X, y = check_X_y(X, y, accept_sparse=("csr", "csc"), multi_output=True,
                         y_numeric=True)

        K = self._get_kernel(X)
        alpha = np.atleast_1d(self.alpha)

        ravel = False
        if len(y.shape) == 1:
            y = y.reshape(-1, 1)
            ravel = True

        copy = self.kernel == "precomputed"
        self.dual_coef_ = _solve_cholesky_kernel(K, y, alpha,
                                                 sample_weight,
                                                 copy)
        if ravel:
            self.dual_coef_ = self.dual_coef_.ravel()

        self.X_fit_ = X

        return self



def _solve_cholesky_kernel(K, y, alpha, sample_weight=None, copy=False):
    # dual_coef = inv(X X^t + alpha*Id) y
    n_samples = K.shape[0]
    n_targets = y.shape[1]

    if copy:
        K = K.copy()

    alpha = np.atleast_1d(alpha)
    one_alpha = (alpha == alpha[0]).all()
    has_sw = isinstance(sample_weight, np.ndarray) \
        or sample_weight not in [1.0, None]

    if has_sw:
        # Unlike other solvers, we need to support sample_weight directly
        # because K might be a pre-computed kernel.
        sw = np.sqrt(np.atleast_1d(sample_weight))
        y = y * sw[:, np.newaxis]
        K *= np.outer(sw, sw)

    if one_alpha:
        # Only one penalty, we can solve multi-target problems in one time.
        K.flat[::n_samples + 1] += alpha[0]

        try:
            # Note: we must use overwrite_a=False in order to be able to
            #       use the fall-back solution below in case a LinAlgError
            #       is raised
            if use_gpu:
                pass

            else:
                dual_coef = linalg.solve(K, y, sym_pos=True,
                                         overwrite_a=False)
            print "Jos is Here!"
        except np.linalg.LinAlgError:
            warnings.warn("Singular matrix in solving dual problem. Using "
                          "least-squares solution instead.")
            dual_coef = linalg.lstsq(K, y)[0]

        # K is expensive to compute and store in memory so change it back in
        # case it was user-given.
        K.flat[::n_samples + 1] -= alpha[0]

        if has_sw:
            dual_coef *= sw[:, np.newaxis]

        return dual_coef
    else:
        # One penalty per target. We need to solve each target separately.
        dual_coefs = np.empty([n_targets, n_samples])

        for dual_coef, target, current_alpha in zip(dual_coefs, y.T, alpha):
            K.flat[::n_samples + 1] += current_alpha

            dual_coef[:] = linalg.solve(K, target, sym_pos=True,
                                        overwrite_a=False).ravel()

            K.flat[::n_samples + 1] -= current_alpha

        if has_sw:
            dual_coefs *= sw[np.newaxis, :]

        return dual_coefs.T

