import unittest
from pyscses.defect_at_site import DefectAtSite
from pyscses.site import Site
from unittest.mock import patch, Mock

class TestDefectAtSite( unittest.TestCase ):

    def setUp( self ):
        self.defect = DefectAtSite( 'A', valence=1.0, mole_fraction=1.0, energy=1.5, site=Mock( spec=Site ), mobility=1 )

    def test_potential_energy( self ):
        phi = 0.5
        self.assertEqual( self.defect.potential_energy( phi ), 2.0 ) 

    def test_boltzmann_one(self):
        phi = 0.5
        temp = 298.0
        with patch('pyscses.defect_at_site.DefectAtSite.potential_energy') as mock_potential_energy:
            mock_potential_energy.return_value = 2.0
            self.assertAlmostEqual(self.defect.boltzmann_one(phi=phi, temp=temp), 1.4996341881113563e-34, places=40)
            mock_potential_energy.assert_called_with(phi)

if __name__ == '__main__':
    unittest.main()
