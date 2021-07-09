import unittest
from unittest.mock import Mock, patch, mock_open, call
from pyscses.set_up_calculation import sites_data_from_file
from pyscses.site_data import SiteData

class TestSetUpCalculation(unittest.TestCase):

    @patch('pyscses.set_up_calculation.SiteData.from_input_string')
    def test_sites_data_from_file_calls_SiteData_from_input_string(self,
        mock_from_input_string):
        sites_input_data = ("A -2.0 1.2345 B -1.0 C 1.0\n"
                            "B +1.0 -0.234 D +0.5\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = 1.2345
        mock_site_data[1].x = -0.234
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat',
                                              x_limits=(-1.0, +1.0))
        expected_calls = [call(line) for line in sites_input_data.split("\n") if line]
        mock_from_input_string.assert_has_calls(expected_calls)

    @patch('pyscses.set_up_calculation.SiteData.from_input_string')
    def test_sites_data_from_file_removes_sites_outside_x_limits(self,
        mock_from_input_string):
        sites_input_data = ("A -2.0 -1.2345 B -1.0 C 1.0\n"
                            "B +1.0 0.234 D +0.5\n"
                            "A -2.0 1.2345  B -1.0 C 1.0\n")
        mock_site_data = [Mock(spec=SiteData) for i in range(3)]
        mock_site_data[0].x = -1.5
        mock_site_data[1].x = 0.0
        mock_site_data[2].x = +1.5
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat',
                                              x_limits=(-1.0, +1.0))
        self.assertEqual(len(sites_data), 1)
        self.assertEqual(sites_data, [mock_site_data[1]])

    @patch('pyscses.set_up_calculation.SiteData.from_input_string')
    def test_sites_data_from_file_returns_SiteData_sorted_by_x(self,
        mock_from_input_string):
        sites_input_data = ("A -2.0 1.2345 B -1.0 C 1.0\n"
                            "B +1.0 -0.234 D +0.5\n")
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].x = +1.0
        mock_site_data[1].x = -1.0
        mock_from_input_string.side_effect = mock_site_data
        with patch('builtins.open', mock_open(read_data=sites_input_data)):
            sites_data = sites_data_from_file(filename='sites.dat',
              x_limits=(-2.0, +2.0))
        self.assertEqual(sites_data, mock_site_data[::-1])


if __name__ == '__main__':
    unittest.main()
