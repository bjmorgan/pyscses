import unittest
from pyscses.set_of_sites import Set_of_Sites

class TestSetOfSites( unittest.TestCase ):

    def test_set_of_sites_is_initialised( self ):
        sites = [ 1,2,3,4 ]
        set_of_sites = Set_of_Sites( sites )
        self.assertEqual( set_of_sites.sites, sites )

    def test_set_of_sites_addition( self ):
        sites1 = [ 1, 2, 3, 4 ]
        sites2 = [ 5, 6, 7, 8 ]
        set_of_sites_1 = Set_of_Sites( sites1 )
        set_of_sites_2 = Set_of_Sites( sites2 )
        combined_set_of_sites = set_of_sites_1 + set_of_sites_2
        self.assertEqual( combined_set_of_sites.sites, sites1 + sites2 )

if __name__ == '__main__':
    unittest.main()
