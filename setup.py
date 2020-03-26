#! /bin/env python

#from distutils.core import setup
from setuptools import setup

setup(name='CINDES',
      version='1.0',
      description='Combinatiorial INverse DESigner',
      author='Jos Teunissen',
      author_email='jteuniss@vub.ac.be',
      scripts=['scripts/cyreader.py', 'scripts/ID_gauss'],  # makes that automatically the right !# path is used.
      #data_files=[ ("", ['requirements.txt'])],
      #py_modules=['foo'],

      # this works only in setuptools not in distutils
      install_requires=[
          "numpy",
          "scipy",
          "sklearn",
          #"pybel",
          #"seaborn"
          ]
     )


# probably want to run as: python setup.py build --build-base=/u/jteuniss/CINDES_setup/test_setup
