''' main callable for INDES

This module is executed when running python -m CINDES [arguments]

Example: nohup python -m CINDES -i inputfile > outputfile 2>&1 &
'''

import argparse
import logging
import sys
import time

from CINDES.utils.writings import print_title
from CINDES.inputreader import read_input

logging.getLogger().setLevel(logging.INFO)


class Unbuffered(object):
    '''to make an unbuffered print interface'''

    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


sys.stdout = Unbuffered(sys.stdout)

if True:
    print time.ctime()
    print_title(
        "C I N D E S\nAn Inverse Molecular Design Program\nwritten by Jos L. Teunissen",
        newlines=True)

    # READ COMMAND LINE ARGUMENTS
    parser = argparse.ArgumentParser(description="INverse DESign package")
    parser.add_argument(
        "-i",
        "--inputfile",
        type=str,
        default='INPUT',
        help="name of the input file. default name: INPUT")
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="increase output verbosity")
    args = parser.parse_args()

    # INPUT READING
    logging.info("name of input-file:" + args.inputfile)
    run = read_input(args.inputfile)

    # START PROGRAM PROCEDURE
    if run.procedure in ['standard', 'bfs']:
        from algorithms.BFS import BFS
        BFS(run)

    elif run.procedure in ['sd', 'rsd'
                           ]:  # can be sd or reduced steepest descent, rsd
        from algorithms.BFS import SD
        SD(run)

    elif run.procedure in ['empty'
                           ]:  # can be sd or reduced steepest descent, rsd
        from algorithms.BFS import runEmpty
        runEmpty(run)

    elif run.procedure in ['ga', 'genetic algorithm', 'genalg']:
        from algorithms import GA
        GA.main(run)

    elif run.procedure in ['deap']:
        from algorithms import DEAP
        DEAP.main(run)

    elif run.procedure in ['pso', 'cpso', 'particleswarm']:
        from algorithms import PSO
        PSO.main(run)

    elif run.procedure in ['gp']:
        from algorithms import GPopt
        GPopt.GP(run)

    elif run.procedure == 'test':
        from algorithms.procedures import testrun
        testrun(run)

    elif run.procedure == 'generate':
        from algorithms.procedures import generate_procedure
        generate_procedure(run)

    elif run.procedure == 'genconf':
        from algorithms.procedures import genconf
        genconf(run)

    elif run.procedure in ['getrandom', 'genrandom']:
        from algorithms.procedures import genrandom
        genrandom(run)

    elif run.procedure == "testpred":
        from algorithms.procedures import testpred
        testpred(run)

    elif run.procedure == 'getdivers':
        from utils.getdivers import database_construction
        database_construction(run)

    else:
        logging.warning('proceduretype not recognized')

    print "EOF __main__.py"
