    def linreg(self,model,alpha=0,intercept=False,printing=0,twosite=False,**kwargs):
        '''does the linear regression and finds the useful parameters'''

        if model == 'LinearRegression':
            clf = linear_model.LinearRegression(fit_intercept=intercept)
        elif model in ['Ridge']:
            clf = linear_model.Ridge(alpha=alpha,fit_intercept=intercept,tol=0.001,solver='auto')
        elif model in ['RidgeCV','ridgecv']:
            clf = linear_model.RidgeCV(alphas=alpha, fit_intercept=intercept, store_cv_values=True)
        elif model in ['Lasso']:
            clf = linear_model.Lasso(alpha=alpha,fit_intercept=intercept,tol=0.001)
        elif model in ['ElasticNet']:
            l1_ratio = 0.1 #default 0.5
            alpha = 1e-2
            clf = linear_model.ElasticNet(alpha=alpha,l1_ratio=l1_ratio, fit_intercept=intercept,tol=0.001)

        if args.verbose>2:
            if twosite:
                sprint(5,self.X2)
            else:
                sprint(5,self.X)
            sprint(5,self.Y)
        if args.fraction:
          if twosite:
            self.X_train, self.X_test, self.Y_train, self.Y_test = train_test_split(self.X2,self.Y, train_size = args.fraction)
            print "size training set:", np.shape(self.Y_train)
            print "size test set:", np.shape(self.Y_test)
            clf.fit(self.X_train, self.Y_train)
          else:
            self.X_train, self.X_test, self.Y_train, self.Y_test = train_test_split(self.X,self.Y, train_size = args.fraction)
            print "size training set:", np.shape(self.Y_train)
            print "size test set:", np.shape(self.Y_test)
            clf.fit(self.X_train, self.Y_train)
        else:
            if twosite:
                clf.fit(self.X2,self.Y)
            else:
                clf.fit(self.X,self.Y)

        return clf

    def predict(self,indices,clf):
        pre_confs = [ indtocon(item) for item in indices ]
        pre_Xs = self.extractX(pre_confs)
        #pre_2Xs= self.extract2D(pre_confs)
        print "prediction. 2 Xs:"
        preds = clf.predict(pre_Xs)
        sprint(2,pre_Xs,preds)
        return preds


    def predict2(self,indices,clf):
        pre_confs = [ indtocon(item) for item in indices ]
        if debug: print "indices:", indices
        #pre_Xs = self.extractX(pre_confs)
        pre_2Xs= self.extract2DX(pre_confs)
        print "two dimensional prediction. 2 Xs:"
        preds = clf.predict(pre_2Xs)
        sprint(2,pre_2Xs,preds)
        return preds


def do():
    for i in range(args.times[0]):
        if args.ols:
            clf_LS = myrun.linreg(model=linmodels[0], intercept=args.intercept)
            if args.intercept:
                print "intercept:", clf_LS.intercept_
            if args.analyse:
                myrun.linreg_analyse(clf_LS)
        if args.ridge:
            print "args.ridge:", args.ridge
            clf_Ridge = myrun.linreg(model='Ridge', alpha=args.ridge, intercept=args.intercept)
            if args.intercept:
                print "intercept:", clf_Ridge.intercept_
            if args.analyse:
                errors = myrun.linreg_analyse(clf_Ridge,model='Ridge')
                allerrors.append(errors)
        if args.ridgecv:
            print "args.ridgeCV"
            alpha = [ 10**i for i in np.arange(-10,10,0.5) ]
            clf_RidgeCV = myrun.linreg(model='RidgeCV',alpha=alpha)
            if args.intercept:
                print "intercept:", clf_RidgeCV.intercept_
            if args.analyse:
                myrun.linreg_analyse(clf_RidgeCV,model='RidgeCV')
        if args.lasso:
            print "args.lasso"
            alpha = args.lasso
            clf_Lasso = myrun.linreg(model='Lasso',alpha=alpha)
            if args.analyse:
                myrun.linreg_analyse(clf_Lasso,model='Lasso')
        if args.twosite:
            hits = myrun.extract2()
            print hits
            if True:
                alpha=50
                print "args.ridge:", alpha
                clf_Ridge2D = myrun.linreg(model='Ridge', alpha=alpha,twosite=True)
                if args.analyse:
                    errors = myrun.linreg_analyse2(clf_Ridge2D,hits=hits,model='Ridge')
                    allerrors.append(errors)
            if False:
                alpha=1e-2
                print "args.Lasso:", alpha
                clf_Lasso2D = myrun.linreg(model='Lasso', alpha=alpha,twosite=True)
                if args.analyse:
                    errors = myrun.linreg_analyse2(clf_Lasso2D,hits=hits,model='Lasso')
                    allerrors.append(errors)
            if False:
                alpha=1e-2
                print "args.ElasticNet:", alpha
                clf_EN2D = myrun.linreg(model='ElasticNet', alpha=alpha,twosite=True)
                if args.analyse:
                    errors = myrun.linreg_analyse2(clf_EN2D,hits=hits,model='ElasticNet')
                    allerrors.append(errors)
    if args.ridge or args.ols or args.ridgecv or args.twosite:
        try:
            for item in allerrors:
                print " ".join(map(str,item))
            print "means RMSE_train/MAE_train/RMSE_test/MAE_test:", np.mean(allerrors,axis=0)
        except (ValueError,TypeError):
            print "error error"
            for item in allerrors:
                print item
            pass

    #myrun.difmodel()
    return myrun


#####################################
#####     END MAIN PROGRAM     ######
#####################################


#####################################
#####   START CALL FROM CINDES   ####
#####################################
#for default args:
class defaults(object):
   def __init__(self):
       self.verbose=1
       self.plot=0
       self.twosite=True
       self.xyplot=0
       self.symmetry=False
       self.column=1
       self.analyse=False
       self.equalsites=False
       self.fraction=None
       self.intersect=False
args = defaults()

@log_io()
def twodim_regression(table, indices,identify,column=2, **kwargs):
    print "In call in fitter.py"
    global args
    args.column=column
    if any(item in identify for item in ['ada', 'adhoma']):
        myrun = Adamantane('ada')
    elif any(item in identify for item in ['dia','dilu','diho','dimi','dima','dilumi']):
        myrun = Diamantane('dia')
    else:
        raise SystemExit('No identify_ identified')
    myrun.extract(table=table)
    hits = myrun.extract2()   #different. 
    print hits
    # DETERMINE ALPHA:
    if True:
        alphas = [ 1*10**i for i in [ -4, -2, -1, 0, 1, 2, 4 ] ]
        clf_RidgeCV = myrun.linreg(model='RidgeCV', twosite=True, alpha=alphas)
        alpha = clf_RidgeCV.alpha_
        print "alpha used:", alpha
    else:
        alpha=1e-4
    print "args.ridge alpha parameter:", alpha
    clf_Ridge2D = myrun.linreg(model='Ridge', alpha=alpha,twosite=True)
    errors = myrun.linreg_analyse2(clf_Ridge2D, model='Ridge')
    predictions = myrun.predict2(indices,clf_Ridge2D)
    print "two_dimensional predictions:", predictions
    return predictions

#####################################
#####    END CALL FROM CINDES   #####
#####################################


