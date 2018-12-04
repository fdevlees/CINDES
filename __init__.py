''' the whole CINDES package '''

__version__ = '0.4'
__author__ = 'J.L. Teunissen'

'''
 Here the rootlogger is configured to display only the message
 At the moment there are still many multiple-line logging event that are not displayed
 correctly when other formatters are used such as:
>formatter = logging.Formatter('%(levelname)8s %(name)s | %(message)s')

 There are also still a lot of print statements in the code that have to be
 converted to logging events

 a more advanced system of loggers could eventually be used by
>logger = logging.getLogger(__name__)
 but for the moment this does not seem necessary.
'''
import sys
import logging
logger = logging.getLogger()
logger.setLevel(logging.WARNING)  # This toggles all the logging
ch = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(message)s')
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)
logger.addHandler(ch)

# exclude scripts/tests/old_modules from "from CINDES import *
__all__ = ['INDES', 'cclib', 'pyevolve', 'PSO', 'utils', 'predictor']
from CINDES import *

# move all module in INDES to top namespace
from INDES import *

