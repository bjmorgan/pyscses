import unittest
from pyscses.site_data import SiteData, InputFormatError
from pyscses.site_data import DefectData
from unittest.mock import patch

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

    @patch('pyscses.site_data.SiteData.input_string_is_valid_syntax')
    def test_from_input_string(self, mock_input_string_is_valid_syntax):
        mock_input_string_is_valid_syntax.return_value = True
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

    @patch('pyscses.site_data.SiteData.input_string_is_valid_syntax')
    def test_from_input_string_raise_InputFormatError_if_input_string_is_invalid(self,
        mock_input_string_is_valid_syntax):
        input_string = "A -2.0 1.2345 B -1.0 C 1.0 X"
        mock_input_string_is_valid_syntax.return_value = False
        with self.assertRaises(InputFormatError):
            site_data = SiteData.from_input_string(input_string,
                                                   validate_input=True)
        mock_input_string_is_valid_syntax.assert_called_with(input_string)

    @patch('pyscses.site_data.SiteData.input_string_is_valid_syntax')
    def test_from_input_string_calls_input_string_is_valid_as_default_behaviour(self,
        mock_input_string_is_valid_syntax):
        input_string = "A -2.0 1.2345 B -1.0 C 1.0 X"
        site_data = SiteData.from_input_string(input_string)
        mock_input_string_is_valid_syntax.assert_called_with(input_string)


class TestSiteDataStaticMethods(unittest.TestCase):

        def test_input_string_is_valid_syntax_returns_True(self):
            input_string = "A -2.0 1.2345 B -1.0 C 1.0"
            self.assertEqual(SiteData.input_string_is_valid_syntax(input_string), True)

        def test_input_string_is_valid_syntax_returns_False(self):
            input_string = "A -2.0 1.2345 B -1.0 C 1.0 X"
            self.assertEqual(SiteData.input_string_is_valid_syntax(input_string), False)

        def test_input_string_is_valid_syntax_returns_true_for_valid_string_with_extended_whitespace(self):
            input_string = "A -2.0 1.2345  B -1.0 C 1.0"
            self.assertEqual(SiteData.input_string_is_valid_syntax(input_string), True)

        def test_input_string_is_valid_syntax_retuns_False_example_2(self):
            input_string = "B +1.0 -0.234 D +0.5 E"
            self.assertEqual(SiteData.input_string_is_valid_syntax(input_string), False)

        def test_input_string_is_valid_syntax_returns_True_example_2(self):
            input_string = "O -2 -1.0958833540000002e-09 Vo 8.73134e-05"
            self.assertTrue(SiteData.input_string_is_valid_syntax(input_string))




if __name__ == '__main__':
    unittest.main()
