import unittest
from pyscses.defect_at_site import DefectAtSite
from pyscses.site import Site
from unittest.mock import patch, Mock

class TestDefectAtSiteInit(unittest.TestCase):
    
    def test_init(self):
        mock_site = Mock(spec=Site)
        defect_at_site = DefectAtSite(label='A',
                                      valence=1.0,
                                      mole_fraction=0.5,
                                      mobility=2.3,
                                      energy=-1.5,
                                      site=mock_site)
        self.assertEqual(defect_at_site.label, 'A')
        self.assertEqual(defect_at_site.valence, 1.0)
        self.assertEqual(defect_at_site.mole_fraction, 0.5)
        self.assertEqual(defect_at_site.mobility, 2.3)
        self.assertEqual(defect_at_site.site, mock_site),
        self.assertEqual(defect_at_site.fixed, False)
        
    def test_init_with_fixed_equals_true(self):
        mock_site = Mock(spec=Site)
        defect_at_site = DefectAtSite(label='A',
                                      valence=1.0,
                                      mole_fraction=0.5,
                                      mobility=2.3,
                                      energy=-1.5,
                                      site=mock_site,
                                      fixed=True)
        self.assertEqual(defect_at_site.fixed, True)
        
class TestDefectAtSite(unittest.TestCase):

    def setUp(self):
        self.defect = DefectAtSite('A',
                                   valence=1.0,
                                   mole_fraction=1.0,
                                   energy=1.5,
                                   site=Mock(spec=Site),
                                   mobility=1)

    def test_potential_energy(self):
        phi = 0.5
        self.assertEqual(self.defect.potential_energy(phi), 2.0) 

    def test_boltzmann_one(self):
        phi = 0.5
        temp = 298.0
        with patch('pyscses.defect_at_site.DefectAtSite.potential_energy') as mock_potential_energy:
            mock_potential_energy.return_value = 2.0
            self.assertAlmostEqual(self.defect.boltzmann_one(phi=phi, temp=temp), 
                                   1.4996341881113563e-34, places=40)
            mock_potential_energy.assert_called_with(phi)
            
    def test_boltzmann_two(self):
        phi = 0.5
        temp = 298.0
        self.defect.mole_fraction = 1.5
        with patch('pyscses.defect_at_site.DefectAtSite.boltzmann_one') as mock_boltzmann_one:
            mock_boltzmann_one.return_value = 2.0
            self.assertEqual(self.defect.boltzmann_two(phi, temp), 3.0)
            mock_boltzmann_one.assert_called_with(phi, temp)
            
    def test_boltzmann_three(self):
        phi = 0.5
        temp = 298.0
        self.defect.mole_fraction = 1.5
        with patch('pyscses.defect_at_site.DefectAtSite.boltzmann_one') as mock_boltzmann_one:
            mock_boltzmann_one.return_value = 2.0
            self.assertEqual(self.defect.boltzmann_three(phi, temp), 1.5)
            mock_boltzmann_one.assert_called_with(phi, temp)
        
            

if __name__ == '__main__':
    unittest.main()
