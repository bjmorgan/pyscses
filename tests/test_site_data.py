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

    def test_as_input_line(self):
        expected_string = "A -2.0 1.2345 B -1.0 C 1.0"
        self.assertEqual(self.site_data.as_input_line(),
                         expected_string)


if __name__ == '__main__':
    unittest.main()
