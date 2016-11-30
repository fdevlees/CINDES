
import unittest
#from .context import INDES
import INDES

class TestINDES1(unittest.TestCase):

    def setUp(self):
        self.indices = ['CNOO_CNHH','CNOO_CNOO_CNOO']

    def tearDown(self):
        self.indices = None

    def test_skipper(self):
        out = INDES.skipper(self.indices,[])
        expected_outcome = [['CNOO_CNHH', 1, 72], ['CNOO_CNOO_CNOO', 1, 129]]
        self.assertEqual( out, expected_outcome )

class TestINDES2(unittest.TestCase):
    def setUp(self):
        self.param = {'nrandom': 8}
        self.array = [ [ 'CNOO' , 'NH2' ] , [ 'HH', 'CC', 'CCl' ] , [ 'Cl', 'F', 'Br' ] ]

    def tearDown(self):
        pass

    def test_genrandom(self):
        self.assertTrue( INDES.genrandom(self.param, self.array) )

if __name__ == "__main__":
    unittest.main()


