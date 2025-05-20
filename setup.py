#! /bin/env python

#from distutils.core import setup
from setuptools import setup, find_packages

setup(name='CINDES',
      version='1.0',
      description='Combinatiorial INverse DESigner',
      authors='Jos Teunissen, Eline Desmedt, David Smets',
      author_email='fdevlees@vub.be',
      scripts=['scripts/cyreader.py', 'scripts/ID_gauss'],  # makes that automatically the right !# path is used 
      packages=find_packages(where='CINDES'),
      package_dir={'': 'CINDES'},
      install_requires=[
          "numpy",
          "scipy",
          "scikit-learn"
          ]
     )

