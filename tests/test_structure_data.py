import unittest
from pyscses.site_data import SiteData
from pyscses.structure_data import StructureData, SplitSitesData
from unittest.mock import Mock, patch
import numpy as np

class TestStructureData(unittest.TestCase):

    @patch('pyscses.structure_data.StructureData.split_sites_data')
    def test_structure_data_is_initialised(self,
                                           mock_split_sites_data):
        sites_data = [Mock(spec=SiteData),
                      Mock(spec=SiteData)]
        x_limits = (0.0, 1.0)
        b = 1.0
        c = 1.0
        system = 'single'
        mock_site_data = [Mock(spec=SiteData),
                          Mock(spec=SiteData)]
        mock_site_data[0].defect_labels = {'A'}
        mock_site_data[1].defect_labels = {'B'}
        mock_site_data[0].x = 0.5
        mock_site_data[1].x = 0.7
        mock_split_sites_return = SplitSitesData(mock_site_data,
                                   (Mock(spec=SiteData),
                                    Mock(spec=SiteData)))
        mock_split_sites_data.return_value = mock_split_sites_return
        structure_data = StructureData(sites_data=sites_data,
                                       x_limits=x_limits,
                                       b=b,
                                       c=c,
                                       system=system)
        self.assertEqual(structure_data.sites_data, mock_split_sites_return[0])
        self.assertEqual(structure_data.adjacent_sites_data, mock_split_sites_return[1])
        self.assertEqual(structure_data.x_limits, x_limits)
        self.assertEqual(structure_data.b, b)
        self.assertEqual(structure_data.c, c)
        np.testing.assert_array_equal(structure_data.site_x_coords, np.array([0.5, 0.7]))
        self.assertEqual(structure_data.defect_labels, {'A', 'B'})
        np.testing.assert_array_equal(structure_data.site_x_coords_by_defect['A'], np.array([0.5]))
        np.testing.assert_array_equal(structure_data.site_x_coords_by_defect['B'], np.array([0.7]))
        self.assertEqual(structure_data.system, system)

    def test_split_sites_data(self):
        mock_sites_data = [Mock(spec=SiteData),
                           Mock(spec=SiteData),
                           Mock(spec=SiteData),
                           Mock(spec=SiteData)]
        for sd, x in zip(mock_sites_data, [0.0, 1.0, 2.0, 3.0]):
            sd.x = x
        x_limits = (0.5, 2.5)
        split_sites_data = StructureData.split_sites_data(sites_data=mock_sites_data,
                         x_limits=x_limits)
        self.assertEqual(split_sites_data.inner_sites_data, mock_sites_data[1:3])
        self.assertEqual(split_sites_data.adjacent_sites_data, tuple(mock_sites_data[0::3]))

    def test_split_sites_data_raises_ValueError_for_fewer_than_three_site_coordinates(self):
        mock_sites_data = [Mock(spec=SiteData),
                           Mock(spec=SiteData),
                           Mock(spec=SiteData)]
        mock_sites_data[0].x = 1.0
        mock_sites_data[1].x = 3.0
        mock_sites_data[2].x = 3.0
        x_limits = (0.5, 2.5)
        with self.assertRaises(ValueError):
            split_sites_data = StructureData.split_sites_data(sites_data=mock_sites_data,
                 x_limits=x_limits)

    def test_split_sites_data_raise_ValueError_if_no_site_adjacent_to_lower_x_limit(self):
        mock_sites_data = [Mock(spec=SiteData),
                           Mock(spec=SiteData),
                           Mock(spec=SiteData)]
        mock_sites_data[0].x = 1.0
        mock_sites_data[1].x = 2.0
        mock_sites_data[2].x = 3.0
        x_limits = (0.5, 2.5)
        with self.assertRaises(ValueError):
            split_sites_data = StructureData.split_sites_data(sites_data=mock_sites_data,
                x_limits=x_limits)

    def test_split_sites_data_raise_ValueError_if_no_site_adjacent_to_upper_x_limit(self):
        mock_sites_data = [Mock(spec=SiteData),
                           Mock(spec=SiteData),
                           Mock(spec=SiteData)]
        mock_sites_data[0].x = 0.0
        mock_sites_data[1].x = 1.0
        mock_sites_data[2].x = 2.0
        x_limits = (0.5, 2.5)
        with self.assertRaises(ValueError):
            split_sites_data = StructureData.split_sites_data(sites_data=mock_sites_data,
                x_limits=x_limits)

    @patch('pyscses.structure_data.sites_data_from_file')
    @patch('pyscses.structure_data.StructureData')
    def test_StructureData_from_file(self,
                                     mock_StructureData,
                                     mock_sites_data_from_file):
        filename = 'foo'
        x_limits = (0.0, 1.0)
        b = 1.0
        c = 1.0
        system = 'single'
        mock_sites_data = [Mock(spec=SiteData)]
        mock_sites_data_from_file.return_value = mock_sites_data
        StructureData.from_file(filename=filename,
                                x_limits=x_limits,
                                b=b,
                                c=c,
                                system=system)
        mock_StructureData.assert_called_with(sites_data=mock_sites_data,
                                              x_limits=x_limits,
                                              b=b,
                                              c=c,
                                              system=system)





if __name__ == '__main__':
    unittest.main()
