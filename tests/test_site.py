import unittest
from pyscses.set_of_sites import SetOfSites
from pyscses.defect_species import DefectSpecies
from pyscses.defect_at_site import DefectAtSite
from pyscses.site_data import SiteData
from pyscses.defect_data import DefectData
from pyscses.site import Site, LabelError
from unittest.mock import Mock, patch
from pyscses.constants import fundamental_charge
import numpy as np

def create_mock_defect_species(n):
    labels = ['a', 'b', 'c', 'd', 'e']
    valence = [-2.0, -1.0, 0.0, 1.0, 2.0]
    mole_fraction = [0.15, 0.25, 0.35, 0.45, 0.55]
    mobility = [0.1, 0.2, 0.3, 0.4, 0.5]
    mock_defect_species = []
    for i in range(n):
        m = Mock(spec=DefectSpecies)
        m.label = labels.pop()
        m.mole_fraction = mole_fraction.pop()
        m.valence = valence.pop()
        m.mobility = mobility.pop()
        m.fixed = False
        mock_defect_species.append(m)
    return mock_defect_species

def create_mock_defects_at_site(n):
    labels = ['A', 'B', 'C', 'D', 'E']
    valence = [-2.0, -1.0, 0.0, 1.0, 2.0]
    mole_fraction = [0.15, 0.25, 0.35, 0.45, 0.55]
    mobility = [0.1, 0.2, 0.3, 0.4, 0.5]
    energies = [-0.1, -0.2, -0.3, -0.4, -0.5]
    mock_defects_at_site = []
    for i in range(n):
        m = Mock(spec=DefectAtSite)
        m.label = labels.pop()
        m.valence = valence.pop()
        m.mole_fraction = mole_fraction.pop()
        m.mobility = mobility.pop()
        m.energy = energies.pop()
        m.fixed = False
        mock_defects_at_site.append(m)
    return mock_defects_at_site

class TestSiteInit(unittest.TestCase):

    def test_site_is_initialised(self):
        mock_defect_species = create_mock_defect_species(2)
        mock_defects_at_site = create_mock_defects_at_site(2)
        with patch('pyscses.site.DefectAtSite', autospec=True) as mock_DefectAtSite:
            mock_DefectAtSite.side_effect = mock_defects_at_site
            site = Site(label='A',
                        x=1.5,
                        defect_species=mock_defect_species,
                        defect_energies=[-0.2, +0.2])
        self.assertEqual(site.label, 'A')
        self.assertEqual(site.x, 1.5)
        self.assertEqual(site.defect_species, mock_defect_species)
        self.assertEqual(site.defect_energies, [-0.2, +0.2])
        np.testing.assert_equal(site.scaling, np.array([1.0, 1.0]))
        self.assertEqual(site.valence, 0.0)
        self.assertEqual(site.saturation_parameter, 1.0)
        self.assertEqual(site.fixed_defects, ())
        self.assertEqual(site.mobile_defects, tuple(mock_defects_at_site))
        self.assertEqual(site.alpha, 1.0)

    def test_site_is_initialised_with_optional_args(self):
        mock_defect_species = create_mock_defect_species(2)
        with patch('pyscses.site.DefectAtSite', autospec=True) as mock_DefectAtSite:
            mock_DefectAtSite.side_effect = create_mock_defects_at_site(2)
            site = Site(label='B',
                        x=1.5,
                        defect_species=mock_defect_species,
                        defect_energies=[-0.2, +0.2],
                        scaling=[0.5, 0.4],
                        valence=-2.0,
                        saturation_parameter=0.1)
        self.assertEqual(site.label, 'B')
        self.assertEqual(site.x, 1.5)
        self.assertEqual(site.defect_species, mock_defect_species)
        self.assertEqual(site.defect_energies, [-0.2, +0.2])
        np.testing.assert_equal(site.scaling, np.array([0.5, 0.4]))
        self.assertEqual(site.valence, -2.0)
        self.assertEqual(site.saturation_parameter, 0.1)
        self.assertEqual(site.alpha, 0.1)

    def test_site_init_with_mixed_mobile_and_fixed_defects(self):
        mock_defect_species = create_mock_defect_species(3)
        mock_defects_at_site = create_mock_defects_at_site(3)
        mock_defects_at_site[0].fixed = False
        mock_defects_at_site[0].mole_fraction = 0.4
        mock_defects_at_site[1].fixed = True
        mock_defects_at_site[1].mole_fraction = 0.3
        mock_defects_at_site[2].fixed = True
        mock_defects_at_site[2].mole_fraction = 0.2
        with patch('pyscses.site.DefectAtSite', autospec=True) as mock_DefectAtSite:
            mock_DefectAtSite.side_effect = mock_defects_at_site
            site = Site(label='C',
            x=1.5,
            defect_species=mock_defect_species,
            defect_energies=[-0.2, +0.2, 0.0])
        self.assertEqual(site.fixed_defects, (mock_defects_at_site[1], mock_defects_at_site[2]))
        self.assertEqual(site.mobile_defects[0], mock_defects_at_site[0])
        self.assertEqual(site.alpha, 0.5)

    def test_site_init_data_check_1(self):
        """Checks that initialising a Site object raises a ValueError if n(defect_species) != n(defect_energies)"""
        mock_defect_species = create_mock_defect_species(1)
        with patch('pyscses.site.DefectAtSite', autospec=True) as mock_DefectAtSite:
            with self.assertRaises(ValueError):
                site = Site(label='A',
                            x=1.5,
                            defect_species=mock_defect_species,
                            defect_energies=[-0.2, +0.2])

    def test_site_init_data_check_2(self):
        """Checks that initialising a Site object raises a ValueError if n(defect_species) != n(scaling) (if passed)"""
        mock_defect_species = create_mock_defect_species(2)
        with patch('pyscses.site.DefectAtSite', autospec=True) as mock_DefectAtSite:
            with self.assertRaises(ValueError):
                site = Site(label='A',
                            x=1.5,
                            defect_species=mock_defect_species,
                            defect_energies=[-0.2, +0.2],
                            scaling=[0.5])

class TestSite(unittest.TestCase):

    def setUp(self):
        mock_defect_species = create_mock_defect_species(2)
        mock_defects_at_site = create_mock_defects_at_site(2)
        with patch('pyscses.site.DefectAtSite', autospec=True) as mock_DefectAtSite:
            mock_DefectAtSite.side_effect = mock_defects_at_site
            self.site = Site(label='A',
                        x=1.5,
                        defect_species=mock_defect_species,
                        defect_energies=[-0.2, +0.2])


    def test_defect_with_label(self):
        self.site.defects[0].label = 'foo'
        self.site.defects[1].label = 'bar'
        self.assertEqual(self.site.defect_with_label('foo'), self.site.defects[0])
        self.assertEqual(self.site.defect_with_label('bar'), self.site.defects[1])

    def test_defect_with_label_2(self):
        """Checks that defect_with_label() raises a LabelError if the argument does not match any of the defect labels for this site."""
        self.site.defects[0].label = 'foo'
        self.site.defects[1].label = 'bar'
        with self.assertRaises(LabelError):
            self.site.defect_with_label('banana')

    def test_energies(self):
        self.site.defects[0].energy = -0.2
        self.site.defects[1].energy = +0.2
        self.assertEqual(self.site.energies(), [-0.2, +0.2])

    def test_probabilities_one(self):
        self.site.defects[0].boltzmann_factor = Mock(return_value=0.1)
        self.site.defects[0].mole_fraction = 0.2
        self.site.defects[0].label = 'A'
        self.site.defects[1].boltzmann_factor = Mock(return_value=0.1)
        self.site.defects[1].mole_fraction = 0.1
        self.site.defects[1].label = 'B'
        exp_A = ((0.2*0.1/(1.0+(0.2*(0.1-1.0)+0.1*(0.1-1.0)))))
        exp_B= ((0.1*0.1/(1.0+(0.2*(0.1-1.0)+0.1*(0.1-1.0)))))
        self.assertEqual(self.site.probabilities(phi=1.0,
                                                 temp=298.0),
                         {'A': exp_A, 'B': exp_B})

    def test_probabilities_two(self):
        self.site.defects[0].boltzmann_factor = Mock(return_value=0.1)
        self.site.defects[0].mole_fraction = 0.2
        self.site.defects[0].label = 'A'
        self.site.defects[0].fixed = True
        self.site.alpha = 0.8
        self.site.fixed_defects = (self.site.defects[0],)
        self.site.defects[1].boltzmann_factor = Mock(return_value=0.1)
        self.site.defects[1].mole_fraction = 0.1
        self.site.defects[1].label = 'B'
        self.site.mobile_defects = (self.site.defects[1],)
        exp_A = 0.2
        exp_B= 0.8*((0.1*0.1/(0.8+(0.1*(0.1-1.0)))))
        self.assertEqual(self.site.probabilities(phi=1.0,
                                                 temp=298.0),
                         {'A': exp_A, 'B': exp_B})

    def test_charge(self):
        self.site.probabilities = Mock(return_value={'E': 0.1, 'D': 0.2})
        self.site.defects[0].valence = 1.0
        self.site.defects[1].valence = 2.0
        self.site.scaling = 0.5
        self.site.valence = 1.0
        expected_value = ((1.0*0.1 + 2.0*0.2)*0.5 + 1.0) * fundamental_charge
        self.assertEqual(self.site.charge(phi=1.0, temp=298.0),
                         expected_value)

    def test_from_site_data(self):
        mock_defect_species = create_mock_defect_species(2)
        mock_defect_species[0].label = 'X'
        mock_defect_species[1].label = 'Y'
        mock_site_data = Mock(spec=SiteData)
        mock_site_data.label = 'A'
        mock_site_data.valence = -2.0
        mock_site_data.x = 1.2345
        mock_defect_data = [Mock(spec=DefectData), Mock(spec=DefectData)]
        mock_defect_data[0].label = 'X'
        mock_defect_data[1].label = 'Y'
        mock_defect_data[0].energy = -1.0
        mock_defect_data[1].energy = +1.0
        mock_site_data.defect_data = mock_defect_data
        with patch('pyscses.site.Site', autospec=True) as mock_Site:
            site = Site.from_site_data(site_data=mock_site_data,
                                       defect_species=mock_defect_species)
        mock_Site.assert_called_with(label=mock_site_data.label,
                                     x=mock_site_data.x,
                                     defect_species=mock_defect_species,
                                     defect_energies=[-1.0, +1.0],
                                     valence=-2.0)

    def test_from_site_data_only_passes_necessary_defect_species(self):
        mock_defect_species = create_mock_defect_species(3)
        mock_defect_species[0].label = 'X'
        mock_defect_species[1].label = 'Y'
        mock_defect_species[2].label = 'Z'
        mock_site_data = Mock(spec=SiteData)
        mock_site_data.label = 'A'
        mock_site_data.valence = -2.0
        mock_site_data.x = 1.2345
        mock_defect_data = [Mock(spec=DefectData), Mock(spec=DefectData)]
        mock_defect_data[0].label = 'X'
        mock_defect_data[1].label = 'Y'
        mock_defect_data[0].energy = -1.0
        mock_defect_data[1].energy = +1.0
        mock_site_data.defect_data = mock_defect_data
        with patch('pyscses.site.Site', autospec=True) as mock_Site:
            site = Site.from_site_data(site_data=mock_site_data,
                                       defect_species=mock_defect_species)
        mock_Site.assert_called_with(label=mock_site_data.label,
                                     x=mock_site_data.x,
                                     defect_species=mock_defect_species[:2],
                                     defect_energies=[-1.0, +1.0],
                                     valence=-2.0)

    def test_from_site_data_raises_ValueError_if_necessary_defect_species_are_not_passed(self):
        mock_defect_species = create_mock_defect_species(3)
        mock_defect_species[0].label = 'X'
        mock_defect_species[1].label = 'W'
        mock_site_data = Mock(spec=SiteData)
        mock_site_data.label = 'A'
        mock_site_data.valence = -2.0
        mock_site_data.x = 1.2345
        mock_defect_data = [Mock(spec=DefectData), Mock(spec=DefectData)]
        mock_defect_data[0].label = 'X'
        mock_defect_data[1].label = 'Y'
        mock_defect_data[0].energy = -1.0
        mock_defect_data[1].energy = +1.0
        mock_site_data.defect_data = mock_defect_data
        with patch('pyscses.site.Site', autospec=True) as mock_Site:
            with self.assertRaises(ValueError):
                site = Site.from_site_data(site_data=mock_site_data,
                                           defect_species=mock_defect_species)

if __name__ == '__main__':
    unittest.main()
