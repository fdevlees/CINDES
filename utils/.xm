#!/bin/env python
import base64
import os
#with open(os.path.expanduser('~/programs/CINDES/utils/.xm264'),'r') as f:
with open(os.path.expanduser('~/CINDES/utils/.xm264'),'r') as f:
    t = f.read()
print base64.b64decode(t)
