import numpy as np
cimport numpy as np

DTYPE = np.int
ctypedef np.int_t DTYPE_t

def cython_sum(np.ndarray[DTYPE_t, ndim=1] y):
    cdef int N = y.shape[0]
    cdef int x = y[0]
    cdef int i
    for i in xrange(1,N):
        x += y[i]
    return x
