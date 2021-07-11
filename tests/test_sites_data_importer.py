import unittest
from unittest.mock import Mock, patch, mock_open, call
from pyscses.site_data import SiteData, InputFormatError
from pyscses.sites_data_importer import (sites_data_from_file,
    cluster_similar_sites_data)


class TestSetUpCalculation(unittest.TestCase):

    @patch('pyscses.sites_data_importer.SiteData.from_input_string')
    @patch('pyscses.sites_data_importer.cluster_similar_sites_data')
    def test_sites_data_from_file_calls_SiteData_from_input_string(self,
        mock_cluster_similar_sites_data,
        mock_from_input_string):
        sites_input_data = ("A -2.0 1.2345 B -1.0 C 1.0\n"
                            "B +1.0 -0.234 D +0.5\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = 1.2345
        mock_site_data[1].x = -0.234
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat')
        expected_calls = [call(line, validate_input=False)
                          for line in sites_input_data.split("\n") if line]
        mock_from_input_string.assert_has_calls(expected_calls)

    @patch('pyscses.sites_data_importer.SiteData.from_input_string')
    @patch('pyscses.sites_data_importer.cluster_similar_sites_data')
    def test_sites_data_from_file_returns_SiteData_sorted_by_x(self,
        mock_cluster_similar_sites_data,
        mock_from_input_string):
        sites_input_data = ("A -2.0 1.2345 B -1.0 C 1.0\n"
                            "B +1.0 -0.234 D +0.5\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = +1.0
        mock_site_data[1].x = -1.0
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat')
        self.assertEqual(sites_data, mock_site_data[::-1])

    @patch('pyscses.sites_data_importer.SiteData.from_input_string')
    @patch('pyscses.sites_data_importer.SiteData.input_string_is_valid_syntax')
    @patch('pyscses.sites_data_importer.cluster_similar_sites_data')
    def test_sites_data_from_file_validates_input_data(self,
        mock_cluster_similar_sites_data,
        mock_input_string_is_valid_syntax,
        mock_from_input_string):
        mock_input_string_is_valid_syntax.return_value = True
        sites_input_data = ("A -2.0 1.2345 B -1.0 C 1.0\n"
                            "B +1.0 -0.234 D +0.5\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = 1.2345
        mock_site_data[1].x = -0.234
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat')
        expected_calls = [call(line) for line in sites_input_data.split("\n") if line]
        mock_input_string_is_valid_syntax.assert_has_calls(expected_calls)

    @patch('pyscses.sites_data_importer.SiteData.from_input_string')
    @patch('pyscses.sites_data_importer.SiteData.input_string_is_valid_syntax')
    @patch('pyscses.sites_data_importer.cluster_similar_sites_data')
    def test_sites_data_from_file_raises_InputFormatError(self,
        mock_cluster_similar_sites_data,
        mock_input_string_is_valid_syntax,
        mock_from_input_string):
        mock_input_string_is_valid_syntax.return_value = False
        sites_input_data = ("A -2.0 1.2345 B -1.0 C 1.0\n"
                            "B +1.0 -0.234 D +0.5 E\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = 1.2345
        mock_site_data[1].x = -0.234
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            with self.assertRaises(InputFormatError):
                sites_data = sites_data_from_file(filename='sites.dat')

    @patch('pyscses.sites_data_importer.SiteData.from_input_string')
    @patch('pyscses.sites_data_importer.SiteData.input_string_is_valid_syntax')
    @patch('pyscses.sites_data_importer.cluster_similar_sites_data')
    def test_sites_data_from_file_clusters_similar_sites(self,
        mock_cluster_similar_sites_data,
        mock_input_string_is_valid_syntax,
        mock_from_input_string):
        mock_input_string_is_valid_syntax.return_value = True
        sites_input_data = ("A -2.0 1.2e-9 B -1.0 C 1.0\n"
                            "B +1.0 1.25e-9 D +0.5 E\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = 1.2e-9
        mock_site_data[1].x = 1.25e-9
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat')
        mock_cluster_similar_sites_data.assert_called_with(sites_data=mock_site_data,
            distance_threshold=1e-10)

    @patch('pyscses.sites_data_importer.SiteData.from_input_string')
    @patch('pyscses.sites_data_importer.SiteData.input_string_is_valid_syntax')
    @patch('pyscses.sites_data_importer.cluster_similar_sites_data')
    def test_sites_data_from_file_uses_explicit_charges_if_site_charge_is_true(self,
        mock_cluster_similar_sites_data,
        mock_input_string_is_valid_syntax,
        mock_from_input_string):
        mock_input_string_is_valid_syntax.return_value = True
        sites_input_data = ("A -2.0 1.2e-9 B -1.0 C 1.0\n"
                            "B +1.0 1.25e-9 D +0.5 E\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = 1.2e-9
        mock_site_data[1].x = 1.25e-9
        mock_site_data[0].valence = -2.0
        mock_site_data[1].valence = +1.0
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat',
                                              site_charge=True)
        self.assertEqual(sites_data[0].valence, -2.0)
        self.assertEqual(sites_data[1].valence, +1.0)

    @patch('pyscses.sites_data_importer.SiteData.from_input_string')
    @patch('pyscses.sites_data_importer.SiteData.input_string_is_valid_syntax')
    @patch('pyscses.sites_data_importer.cluster_similar_sites_data')
    def test_sites_data_from_file_sets_site_charges_to_zero_if_site_charge_is_false(self,
        mock_cluster_similar_sites_data,
        mock_input_string_is_valid_syntax,
        mock_from_input_string):
        mock_input_string_is_valid_syntax.return_value = True
        sites_input_data = ("A -2.0 1.2e-9 B -1.0 C 1.0\n"
                            "B +1.0 1.25e-9 D +0.5 E\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = 1.2e-9
        mock_site_data[1].x = 1.25e-9
        mock_site_data[0].valence = -2.0
        mock_site_data[1].valence = +1.0
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat',
                                              site_charge=False)
        self.assertEqual(sites_data[0].valence, 0.0)
        self.assertEqual(sites_data[1].valence, 0.0)

    @patch('pyscses.sites_data_importer.SiteData.from_input_string')
    @patch('pyscses.sites_data_importer.SiteData.input_string_is_valid_syntax')
    @patch('pyscses.sites_data_importer.cluster_similar_sites_data')
    def test_sites_data_from_file_sets_site_charges_to_zero_as_default(self,
        mock_cluster_similar_sites_data,
        mock_input_string_is_valid_syntax,
        mock_from_input_string):
        mock_input_string_is_valid_syntax.return_value = True
        sites_input_data = ("A -2.0 1.2e-9 B -1.0 C 1.0\n"
                            "B +1.0 1.25e-9 D +0.5 E\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = 1.2e-9
        mock_site_data[1].x = 1.25e-9
        mock_site_data[0].valence = -2.0
        mock_site_data[1].valence = +1.0
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat')
        self.assertEqual(sites_data[0].valence, 0.0)
        self.assertEqual(sites_data[1].valence, 0.0)

    def test_cluster_similar_sites_data(self):
        mock_sites_data = [Mock(spec=SiteData),
                           Mock(spec=SiteData),
                           Mock(spec=SiteData),
                           Mock(spec=SiteData)]
        mock_sites_data[0].x = 1.0e-10
        mock_sites_data[1].x = 1.9e-10
        mock_sites_data[2].x = 3.0e-10
        mock_sites_data[3].x = 3.9e-10
        cluster_similar_sites_data(mock_sites_data)
        self.assertEqual(mock_sites_data[0].x, 1.45e-10)
        self.assertEqual(mock_sites_data[1].x, 1.45e-10)
        self.assertEqual(mock_sites_data[2].x, 3.45e-10)
        self.assertEqual(mock_sites_data[3].x, 3.45e-10)

    def test_cluster_similar_sites_data_preserves_site_order(self):
        mock_sites_data = [Mock(spec=SiteData),
                           Mock(spec=SiteData),
                           Mock(spec=SiteData),
                           Mock(spec=SiteData)]
        mock_sites_data[0].x = 1.0e-10
        mock_sites_data[1].x = 3.9e-10
        mock_sites_data[2].x = 3.0e-10
        mock_sites_data[3].x = 1.9e-10
        cluster_similar_sites_data(mock_sites_data)
        self.assertEqual(mock_sites_data[0].x, 1.45e-10)
        self.assertEqual(mock_sites_data[1].x, 3.45e-10)
        self.assertEqual(mock_sites_data[2].x, 3.45e-10)
        self.assertEqual(mock_sites_data[3].x, 1.45e-10)

    def test_cluster_similar_sites_data_behaves_sensible_when_passed_a_single_site(self):
        mock_sites_data = [Mock(spec=SiteData)]
        mock_sites_data[0].x = 1.0
        cluster_similar_sites_data(mock_sites_data)


if __name__ == '__main__':
    unittest.main()
