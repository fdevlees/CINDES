''' the whole CINDES package '''

__version__ = '0.4'
__author__ = 'J.L. Teunissen'


# exclude scripts/tests/old_modules from "from CINDES import *
__all__ = ['INDES', 'cclib', 'pyevolve', 'PSO', 'utils', 'predictor']
from CINDES import *

# move all module in INDES to top namespace
from INDES import *

# set a logger
import sys
import logging
logger = logging.getLogger(__name__)

ch = logging.StreamHandler(sys.stdout)
#formatter = logging.Formatter('%(levelname)8s %(name)s | %(message)s')
formatter = logging.Formatter('%(message)s')
ch.setFormatter(formatter)
ch.setLevel(logging.INFO)

logger.addHandler(ch)
logger.setLevel(logging.INFO)  # This toggles all the logging
