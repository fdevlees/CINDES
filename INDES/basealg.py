''' some new not yet used stuff for making the algorithm more object oriented '''

class Algorithm(object):

    def __init__(self, param, array):
        self.param     = param
        self.array     = array
        self.bcok      = 0
        self.startconf = get_startconf(self.param, self.array)
        self.myrun     = Run(**param)
        self.table     = set_table( self.myrun)
        self.maximum   = set_maximum(self.myrun, self.table)
        print(self.myrun) #this should print all the class elements via the __str__ function

    def run():
            # ------------------------------------- # 
            # --- HERE THE MAIN LOOP STARTS --- --- #
            # ------------------------------------- #
    count = 1 # so we start counting at 1!
    while True:
        print_title("COUNT: " + str(count),outline='l',signator="-")

        ### set site order in sequence INPUT: param, count
        sequence = get_sequence(count, myrun)

        ### for each site in sequence:
        for l in range(len(sequence)):
            k = sequence[l]
            print_title("k(site)= " + str(k) + " l(nsite)= "+ str(l),outline='l',signator='=')
            if not l == 0 or count > 1: #define new startconfiguration if not first cycle
                # define new starting geometry
                print "maxsite[0]",maxsite[0]
                del startconf
                startconf = zcon.indtocon(maxsite[0])
                print "newconf: ", startconf


            # STEP 1: INDEXMAKER
            #get indices_all and the indices that still need to be calculated
            # if table is correctly formatted all second element item[1]==1. meaning they are ab-initio calculated
            #indices_todo,data_nodo,configurations,indices_all = zcon.indexmaker2(startconf,array,k,table )
            mols_todo, mols_nodo = zcon.classmaker2(startconf,array,k,table, myrun )
            if 1 in myrun.restrictions:
                mols_todo, mols_nodo = restriction1(mols_todo, mols_nodo, myrun )
            print "----- END random start configurations -----"
            print "indices_todo:",mols_todo
            print "data_nodo:", mols_nodo #all item[1]==1 in data_nodo 

            # STEP 2: PREDICTOR
            # perform prescreaning in a predictions. 
            mols_nocal, mols_tocal, made_pred = predictor(myrun, table, mols_todo,mols_nodo, count, array=array, nsite=l)



            # STEP 3: SUBMITTING PART
            if not myrun.nosub==1:
                mols_all = submittingprocedure(mols_tocal,
                                               mols_nocal,
                                               myrun,
                                             **myrun.TZmat     ) # here call submitting procedure
            else: mols_all = skipper(mols_tocal,mols_nocal)
            print "mols_all:",mols_all

            # STEP 4: SORT
            # sort data in same order as allindices:
            # not necessary anymore in molsclass
            #mols_all = sorted(mols_all, key=lambda x:x.Pvalue)

            # STEP 5: UPDATE DATABASE and LOG results of microiteration
            # logs new elements in data to table and tablebin and whole data to cyclesinfo
            table = loggings(mols_all,table,count,k,l, made_pred, tablename = myrun.tablename)

            # STEP 6: UPDATE OPTIMUM STRUCTURE
            # decide what the maximum site is and if the bc if fullfilled
            print "BCOK:", bcok
            maxsite, bcok = testmax(myrun, mols_all, bcok)

            print("--- %s seconds ---" % (time.time() - myrun.starttime))
            print(myrun.currenttime())
        # END LOOP OVER SITES

        #get maximum and test convergence
        maximum, maxsite,converged = runtest(myrun, maximum, maxsite, count, bcok, mctable=table, array = array)
        if converged==1: break
        count +=1
        if count > param['maxiter']:
            print "maxiterations is reached"
            print "maximum is: ", maximum
            break
    # ---------------------------- #
    # ------ END OF LOOPING ------ # 
    # ---------------------------- #
    print "DONE"
    return
