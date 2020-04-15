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

"""Utilities often used by cclib parsers"""


def convertor(value, fromunits, tounits):
    """Convert from one set of units to another.

    Sources:
        NIST 2010 CODATA (http://physics.nist.gov/cuu/Constants/index.html)
        Documentation of GAMESS-US or other programs as noted

    >>> print "%.1f" % convertor(8, "eV", "cm-1")
    64524.8
    """

    _convertor = {
        "Angstrom_to_bohr": lambda x: x * 1.8897261245,
        "bohr_to_Angstrom": lambda x: x * 0.5291772109,
        "cm-1_to_eV":       lambda x: x / 8065.54429,
        "cm-1_to_hartree":  lambda x: x / 219474.6313708,
        "cm-1_to_kcal":     lambda x: x / 349.7550112,
        "cm-1_to_kJmol-1":  lambda x: x / 83.5934722814,
        "cm-1_to_nm":       lambda x: 1e7 / x,
        "eV_to_cm-1":       lambda x: x * 8065.54429,
        "eV_to_hartree":    lambda x: x / 27.21138505,
        "eV_to_kcal":       lambda x: x * 23.060548867,
        "eV_to_kJmol-1":    lambda x: x * 96.4853364596,

        "hartree_to_cm-1":  lambda x: x * 219474.6313708,
        "hartree_to_eV":    lambda x: x * 27.21138505,
        "hartree_to_kcal":  lambda x: x * 627.50947414,
        "hartree_to_kJmol-1":lambda x: x * 2625.4996398,

        "kcal_to_cm-1":     lambda x: x * 349.7550112,
        "kcal_to_eV":       lambda x: x / 23.060548867,
        "kcal_to_hartree":  lambda x: x / 627.50947414,
        "kcal_to_kJmol-1":  lambda x: x * 4.184,

        "kJmol-1_to_cm-1":  lambda x: x * 83.5934722814,
        "kJmol-1_to_eV":    lambda x: x / 96.4853364596,
        "kJmol-1_to_hartree": lambda x: x / 2625.49963978,
        "kJmol-1_to_kcal":  lambda x: x / 4.184,
        "nm_to_cm-1":       lambda x: 1e7 / x,

        # Taken from GAMESS docs, "Further information",
        # "Molecular Properties and Conversion Factors"
        "Debye^2/amu-Angstrom^2_to_km/mol": lambda x: x * 42.255,

        # Conversion for charges and multipole moments.
        "e_to_coulomb":         lambda x: x * 1.602176565 * 1e-19,
        "e_to_statcoulomb":     lambda x: x * 4.80320425  * 1e-10,
        "coulomb_to_e":         lambda x: x * 0.6241509343 * 1e19,
        "statcoulomb_to_e":     lambda x: x * 0.2081943527 * 1e10,
        "ebohr_to_Debye":       lambda x: x * 2.5417462300,
        "ebohr2_to_Buckingham": lambda x: x * 1.3450341749,
        "ebohr2_to_Debye.ang":  lambda x: x * 1.3450341749,
        "ebohr3_to_Debye.ang2": lambda x: x * 0.7117614302,
        "ebohr4_to_Debye.ang3": lambda x: x * 0.3766479268,
        "ebohr5_to_Debye.ang4": lambda x: x * 0.1993134985,

    }

    return _convertor["%s_to_%s" % (fromunits, tounits)] (value)


class PeriodicTable(object):
    """Allows conversion between element name and atomic no.

    >>> t = PeriodicTable()
    >>> t.element[6]
    'C'
    >>> t.number['C']
    6
    >>> t.element[44]
    'Ru'
    >>> t.number['Au']
    79
    """
    def __init__(self):
        self.element = [None,
            'H', 'He',
            'Li', 'Be',
            'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg',
            'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
            'K', 'Ca',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr',
            'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
            'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
            'Cs', 'Ba',
            'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
            'Fr', 'Ra',
            'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
            'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Uub']
        self.number = {}
        for i in range(1, len(self.element)):
            self.number[self.element[i]] = i

def run_once(f):
    def wrapper(*args, **kwargs):
        if not wrapper.has_run:
            wrapper.has_run = True
            return f(*args, **kwargs)
    wrapper.has_run = False
    return wrapper

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class prettyfloat(float):
    def __repr__(self):
        return "%-0.4f" % self

def round_sig( x, sig=8):
    from math import log10, floor
    try:
        return round(x, sig-int(floor(log10(abs(x))))-1)
    except ValueError:
        if not x==0.0: print("ValueError:", x)
        return x

# this function to redirect the output of the optga keyword
import sys
from contextlib import contextmanager
@contextmanager
def custom_redirection(fileobj):
    old = sys.stdout
    sys.stdout = fileobj
    try:
        yield fileobj
    finally:
        sys.stdout = old

def rm_duplicates(seq, nsig=8):
    seen = set()
    seen_add = seen.add
    new_seq = []
    if True: # try to round to numerical precision errors:
        print("    the sequence(scfenergies?) is rounded to max 10 significant digits")
        seq = [ round_sig( item, sig=10 ) for item in seq ]
    for x in seq:
        if x in seen:
            print("    multiples found in sequence. name probably scfenergies! | value: ", x)
        else:
            new_seq.append(x)
            seen_add(x)
    return new_seq

def import_program(program_name):
    if program_name=='gaussian':
        import gaussian as program
    elif program_name=='orca':
        import orca as program
    elif program_name=='nwchem':
        import nwchem as program
    else:
        print("program_name:", program_name)
        raise SystemExit('not implemented')
    return program

import os
import sys
import traceback
from functools import wraps
from multiprocessing import Process, Queue

def processify(func):
    '''Decorator to run a function as a process.
    Be sure that every argument and the return value
    is *pickable*.
    The created process is joined, so the code does not
    run in parallel.
    '''

    def process_func(q, *args, **kwargs):
        try:
            ret = func(*args, **kwargs)
        except Exception:
            ex_type, ex_value, tb = sys.exc_info()
            error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
            ret = None
        else:
            error = None

        q.put((ret, error))

    # register original function with different name
    # in sys.modules so it is pickable
    process_func.__name__ = func.__name__ + 'processify_func'
    setattr(sys.modules[__name__], process_func.__name__, process_func)

    @wraps(func)
    def wrapper(*args, **kwargs):
        q = Queue()
        p = Process(target=process_func, args=[q] + list(args), kwargs=kwargs)
        p.start()
        ret, error = q.get()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (ex_value.message, tb_str)
            raise ex_type(message)

        return ret
    return wrapper

import string
class PartialFormatter(string.Formatter):
    def __init__(self, missing='~~', bad_fmt='!!'):
        self.missing, self.bad_fmt=missing, bad_fmt

    def get_field(self, field_name, args, kwargs):
        # Handle a key not found
        try:
            val=super(PartialFormatter, self).get_field(field_name, args, kwargs)
            # Python 3, 'super().get_field(field_name, args, kwargs)' works
        except (KeyError, AttributeError):
            val=None,field_name
        return val

    def format_field(self, value, spec):
        # handle an invalid format
        if value==None: return self.missing
        try:
            return super(PartialFormatter, self).format_field(value, spec)
        except ValueError:
            if self.bad_fmt is not None: return self.bad_fmt
            else: raise

def jsonify(data):
    json_data = dict()
    for key, value in data.items():
        if isinstance(value, list):
            value = [ jsonify(item) if isinstance(item, dict) else item for item in value ]
        if isinstance(value, dict):
            value = jsonify(value)
        if isinstance(key, int):
            key = str(key)
        if type(value).__module__=='numpy':
            value = value.tolist()
        json_data[key] = value
    return json_data

def pythonify(json_data):
    for key, value in json_data.items():
        if isinstance(value, list):
            value = [ pythonify(item) if isinstance(item, dict) else item for item in value ]
        elif isinstance(value, dict):
            value = pythonify(value)
        try:
            newkey = int(key)
            del json_data[key]
            key = newkey
        except TypeError:
            pass
        json_data[key] = value
    return json_data


# custom function evaluation
def skipper(mols_tocal, mols_nocal, myrun, iprint=True):
    ''' generate random data '''
    mols_all = mols_tocal + mols_nocal
    for molecule in mols_all:
        if myrun.property=='func':
            # evaluate function and set to arguments
            Pvalue = myrun.function(molecule)
            result = { func_arg:value for func_arg,value in zip(myrun.func_args, Pvalue) }
            print('result', result, Pvalue)
            molecule.props.update(result)
            if isinstance(Pvalue, int) or isinstance(Pvalue, float):
                molecule.Pvalue = Pvalue
            elif len(Pvalue)==1:
                molecule.Pvalue = Pvalue[0]
            else:
                molecule.Pvalue = Pvalue
            #molecule.Pvalue = molecule.props[myrun.func_args[0]]
        else:
            item = molecule.index
            output = 0
            replaced = item.replace('_', '')
            replaced = [x for x in replaced if x.isalpha()]
            for i in replaced:
                try:
                    output += string.uppercase.index(i)
                except ValueError:
                    output += string.lowercase.index(i)
            propx = float(output)
            try:
                if 'bcprop' in param:
                    propy = len(item.replace('_', ''))
                    molecule.boundaries = [float(propy)]
            except NameError:
                pass
            molecule.Pvalue = propx
            molecule.props[myrun.property] = propx
        molecule.predicted = False

    mols_all = mols_tocal + mols_nocal
    return mols_all

class SessionID(object):
    def __init__(self, ID):
        if isinstance(ID, int):
            n_calc = int(str(ID)[1:])
            n_session = int(str(ID)[0])
            name = "{:d}.{:d}".format(n_session, n_calc)
        else:
            name = ID
        self.names= [name]
        self.names.append(float(name))
        splitted = name.split('.')
        zerod = '{}0{}'.format(*splitted)
        numbered = int(zerod)
        self.names.extend([splitted, list(map(int, splitted)), zerod, numbered])
        return

    def __eq__(self, other):
        return other in self.names

    def __repr__(self):
        return "<SessionID object: {}>".format(self.names[0])


if __name__ == "__main__":
    import doctest, utils
    doctest.testmod(utils, verbose=False)
