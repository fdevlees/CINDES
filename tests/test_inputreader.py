''' this test the reader.py module '''
debug=True

import unittest
from CINDES import inputreader as inr
path = '/u/jteuniss/CINDES/tests'

test_param1 = {'cutoff': 0, 'nosub': 2, 'timestep': 300, 'functional': 'b3lyp', 'no1sub': 0, 'tdregression': 0, 'identify': 'neop_', 'semiempirical': 0, 'charge': 0, 'bcprop': 'homo', 'montecarlo': 0, 'program': 'gaussian', 'bcval': '-3', 'ncore': 4, 'maxiter': 10, 'restart': 0, 'test_ready': 2, 'ml': 0, 'timelimit': 25000000, 'multiplejobs': 0, 'property': 'energy', 'procedure': 'standard', 'difmodel': 0, 'nprocs': 2, 'sequence': [], 'ip': 0, 'ea': 0, 'nch3': 4, 'optimum': 'minimum', 'regression': 0, 'extrawaittime': 2, 'line1': [5, 9], 'nosub_file:': '', 'try_ready': 0, 'polar': 0, 'twojob': 0, 'debug': False, 'nrandsites': 2, 'stab': 0, 'basisset': '6-31g', 'mult': 1, 'symlinks': [], 'norandom': 0, 'bcoptimum': 'min', 'startind': ''}

test_param2 = {'cutoff': 0, 'nosub': 2, 'timestep': 300, 'functional': 'b3pw91', 'no1sub': 0, 'tdregression': 0, 'identify': 'neop_', 'semiempirical': 0, 'charge': 0, 'bcprop': 'polar', 'montecarlo': 0, 'program': 'orca', 'bcval': '60', 'ncore': 4, 'maxiter': 10, 'restart': 0, 'test_ready': 1, 'ml': 0, 'timelimit': 25000000, 'multiplejobs': 0, 'property': 'stab', 'procedure': 'standard', 'difmodel': 0, 'nprocs': 5, 'sequence': [0, 1], 'ip': 0, 'sequences': [[0, 1], [1, 0], [1, 0]], 'ea': 0, 'nch3': 4, 'optimum': 'max', 'regression': 1, 'extrawaittime': 2, 'line1': [5, 9, 17], 'nosub_file:': '', 'try_ready': 1, 'polar': 0, 'twojob': 0, 'debug': True, 'nrandsites': 2, 'stab': 1, 'basisset': '6-31g', 'mult': 2, 'symlinks': [], 'norandom': 0, 'bcoptimum': 'min', 'startind': 'N_CF_CNHCHHH'}

test_array2 = [[['N'], ['C', 'F'], ['C', 'Cl']],
               [['C', 'Ph'], ['C', 'F'], ['C', 'Cl']],
                [['C', 'C', 'F', 'H', 'F'],
                ['C', 'N', 'O', 'O'], ['C', 'N', 'H', 'C', 'H', 'H', 'H'], ['C', 'C', 'N']]]

class Test_inputreader1(unittest.TestCase):
    def setUp(self):
        self.inputfile = path + '/testfiles/' + 'INPUT1'
        self.fid = open( self.inputfile, 'r' )
        # to make the printing more verbose when there is an error:
        self.maxDiff = None

    def tearDown(self):
        self.fid.close()

    def test_readfile1(self):
        param1 = inr.readfile(self.fid)
        self.assertEqual( param1, test_param1)

class Test_inputreader2(unittest.TestCase):
    def setUp(self):
        self.inputfile = path + '/testfiles/' + 'INPUT2'
        self.fid = open( self.inputfile, 'r' )
        # to make the printing more verbose when there is an error:
        self.maxDiff = None

    def tearDown(self):
        self.fid.close()

    def test_readfile1(self):
        param2 = inr.readfile(self.fid)
        self.assertEqual( param2, test_param2)

class Test_substireader1(unittest.TestCase):
    def setUp(self):
        self.inputfile = path + '/testfiles/' + 'INPUT2'
        self.fid = open( self.inputfile, 'r' )
        inr.readfile(self.fid)
        # to make the printing more verbose when there is an error:
        self.maxDiff = None

    def tearDown(self):
        self.fid.close()

    def test_substireader(self):
        array2 = inr.substireader( 4, self.fid)
        self.assertEqual( array2, test_array2)

if __name__ == "__main__":
    unittest.main()


