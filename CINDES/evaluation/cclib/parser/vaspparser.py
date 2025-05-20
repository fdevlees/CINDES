# -*- coding: utf-8 -*-
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Parser for VASP output files"""



import re

import numpy

from . import logfileparser
from . import utils


class VASP(logfileparser.Logfile):
    """A VASP ??? file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(VASP, self).__init__(logname="VASP", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "VASP ??? file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'VASP("%s")' % (self.filename)

    def before_parsing(self):
        pass

    def after_parsing(self):
        pass

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""
        pass


if __name__ == "__main__":
    import doctest
    from . import vaspparser
    import sys

    if len(sys.argv) == 1:
        doctest.testmod(vaspparser, verbose=False)

    if len(sys.argv) >= 2:
        parser = vaspparser.VASP(sys.argv[1])
        data = parser.parse()

    if len(sys.argv) > 2:
        for i in range(len(sys.argv[2:])):
            if hasattr(data, sys.argv[2 + i]):
                print(getattr(data, sys.argv[2 + i]))

