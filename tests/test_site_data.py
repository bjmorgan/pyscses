import unittest
from pyscses.site_data import SiteData
from pyscses.site_data import DefectData

class TestDefectDataInit(unittest.TestCase):

    def test_DefectData_is_initialised(self):
        label = 'X'
        energy = -1.0
        defect_data = DefectData(label=label,
                                 energy=energy)
        self.assertEqual(defect_data.label, label)
        self.assertEqual(defect_data.energy, energy)

class TestDefectData(unittest.TestCase):

    def test_eq_returns_true(self):
        label = 'X'
        energy = -1.0
        defect_data_1 = DefectData(label=label,
                                   energy=energy)
        defect_data_2 = DefectData(label=label,
                                   energy=energy)
        self.assertTrue(defect_data_1 == defect_data_2)

    def test_eq_returns_false_on_unequal_labels(self):
        label1 = 'X'
        energy = -1.0
        label2 = 'Y'
        defect_data_1 = DefectData(label=label1,
                                   energy=energy)
        defect_data_2 = DefectData(label=label2,
                                   energy=energy)
        self.assertFalse(defect_data_1 == defect_data_2)

    def test_eq_returns_false_on_unequal_energies(self):
        label = 'X'
        energy1 = -1.0
        energy2 = +1.0
        defect_data_1 = DefectData(label=label,
                                   energy=energy1)
        defect_data_2 = DefectData(label=label,
                                   energy=energy2)
        self.assertFalse(defect_data_1 == defect_data_2)

    def test_eq_evalutates_as_False_for_invalid_comparison(self):
        label = 'X'
        energy = -1.0
        defect_data = DefectData(label=label,
                                 energy=energy)
        self.assertFalse(defect_data == 'foo')


class TestSiteDataInit(unittest.TestCase):

    def test_SiteData_is_initialised(self):
        label = 'A'
        valence = -2.0
        x = 1.2345
        defect_data = (DefectData(label='B',
                                  energy=-1.0),
                       DefectData(label='C',
                                  energy=+1.0))
        site_data = SiteData(label=label,
                             valence=valence,
                             x=x,
                             defect_data=defect_data)
        self.assertEqual(site_data.label, label)
        self.assertEqual(site_data.valence, valence)
        self.assertEqual(site_data.x, x)
        self.assertEqual(site_data.defect_data, defect_data)

class TestSiteData(unittest.TestCase):

    def setUp(self):
        label = 'A'
        valence = -2.0
        x = 1.2345
        defect_data = (DefectData(label='B',
                                  energy=-1.0),
                       DefectData(label='C',
                                  energy=+1.0))
        self.site_data = SiteData(label=label,
                                  valence=valence,
                                  x=x,
                                  defect_data=defect_data)

    def test_as_input_string(self):
        expected_string = "A -2.0 1.2345 B -1.0 C 1.0"
        self.assertEqual(self.site_data.as_input_string(),
                         expected_string)

    def test_from_input_string(self):
        input_string = "A -2.0 1.2345 B -1.0 C 1.0"
        site_data = SiteData.from_input_string(input_string)
        self.assertEqual(site_data.label, 'A')
        self.assertEqual(site_data.valence, -2.0)
        self.assertEqual(site_data.x, 1.2345)
        expected_defect_data = (DefectData(label='B',
                                           energy=-1.0),
                                DefectData(label='C',
                                           energy=+1.0))
        self.assertEqual(site_data.defect_data,
                         expected_defect_data)
                         
class TestSiteDataStaticMethods(unittest.TestCase):
    
        def test_validate_input_string_returns_True(self):
            input_string = "A -2.0 1.2345 B -1.0 C 1.0"
            self.assertTrue(SiteData.validate_input_string(input_string))
            
        def test_validate_input_string_returns_False(self):
            input_string = "A -2.0 1.2345 B -1.0 C 1.0 X"
            self.assertFalse(SiteData.validate_input_string(input_string))
    


if __name__ == '__main__':
    unittest.main()
