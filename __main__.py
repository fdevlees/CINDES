''' main callable for INDES '''
# call with python -m CINDES4 -i <inputfile> 
# nohup python -m CINDES4 -i inputfile > outputfile 2>&1 &

import argparse
import logging
import sys

from utils.writings import print_title
from INDES.inputreader import read_input


class Unbuffered(object):
    def __init__(self,stream):
        self.stream = stream
    def write(self,data):
        self.stream.write(data)
        self.stream.flush()
    def __getattr__(self,attr):
        return getattr(self.stream, attr)
sys.stdout = Unbuffered(sys.stdout)

if True:
    print_title("C I N D E S\nAn Inverse Molecular Design Program\nwritten by Jos L. Teunissen", newlines=True)

    # READ COMMAND LINE ARGUMENTS
    parser = argparse.ArgumentParser(description="INverse DESign package")
    parser.add_argument("-i","--inputfile",type = str,default='INPUTBC',help="name of the input file. default name: INPUTBC")
    parser.add_argument("-z","--zmatrixfile",type = str,default='ZMAT',help="name of the zmatrix file. default name: ZMAT")
    parser.add_argument("-v","--verbose", action="count", default=0, help="increase output verbosity")
    args=parser.parse_args()
    #zmatrixfile is a global variable
    logging.info("name of zmatfile:  " + args.zmatrixfile)
    logging.info("name of input-file:" + args.inputfile)
    # INPUT READING
    param, array = read_input(args.inputfile)
    param['zmatrixfile']=args.zmatrixfile
    # END INPUT READING

    #START PROGRAM PROCEDURE
    if param['procedure'] == 'standard':
        from INDES.procedures import BFS
        BFS(param,array)

    elif 'steepest' in param['procedure']: # can be steepest1 or steepest2
        from INDES.procedures import SteepestDescent
        SteepestDescent(param,array)

    elif param['procedure'] in [ 'ga' , 'genetic algorithm', 'genalg' ]:
        from INDES import GA
        GA.main(param,array)

    elif param['procedure'] == 'test':
        from INDES.procedures import testrun
        testrun(param,array)

    elif param['procedure'] == 'generate':
        from INDES.procedures import generate_procedure
        generate_procedure(param,array)

    elif param['procedure'] == 'genconf':
        from INDES.procedures import genconf
        genconf(param)

    elif param['procedure'] in [ 'getrandom' ,'genrandom']:
        from INDES.procedures import genrandom
        genrandom(param,array)

    else:
        logging.warning('proceduretype not recognized')

    print "EOF __main__.py"
