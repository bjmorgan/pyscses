import unittest
from unittest.mock import Mock, patch, mock_open, call
from pyscses.set_up_calculation import sites_from_file
from pyscses.site_data import SiteData

class TestSetUpCalculation(unittest.TestCase):

    def test_sites_from_file_calls_SiteData_from_input_string(self):
        sites_data = ("A -2.0 1.2345 B -1.0 C 1.0\n"
                      "B +1.0 -0.234 D +0.5\n")
        with patch('builtins.open', mock_open(read_data=sites_data)):
            mock_open.return_value = 'foo'
            with patch('pyscses.set_up_calculation.SiteData.from_input_string') as mock_from_input_string:
                mock_site_data = [Mock(spec=SiteData), Mock(spec=SiteData)]
                mock_from_input_string.return_value = mock_site_data
                sites = sites_from_file(filename='sites.dat',
                                        x_limits=(-1.0, +1.0))
        expected_calls = [call(line) for line in sites_data.split("\n") if line]
        mock_from_input_string.assert_has_calls(expected_calls)


if __name__ == '__main__':
    unittest.main()
