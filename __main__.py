''' main callable for INDES '''

import argparse
import logging

from utils.writings import print_title
from INDES.CINDES3 import read_input, main, testrun, generate_procedure, genconf, genrandom

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
        main(param,array)
    elif param['procedure'] == 'test':
        testrun(param,array)
    elif param['procedure'] == 'generate':
        generate_procedure(param,array)
    elif param['procedure'] == 'genconf':
        genconf(param)
    elif param['procedure'] in [ 'getrandom' ,'genrandom']:
        genrandom(param,array)
    else:
        logging.warning('proceduretype not recognized')
    print "bla"
