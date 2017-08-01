import unittest
from project.grid import Grid, Grid_Point, volumes_from_grid, delta_x_from_grid
from project.set_of_sites import Set_of_Sites
from project.site import Site
from unittest.mock import Mock, MagicMock, patch, call
import numpy as np

class TestGrid( unittest.TestCase ):

    @patch( 'project.grid.volumes_from_grid' )
    @patch( 'project.grid.index_of_grid_at_x' )
    @patch( 'project.grid.Grid_Point' )
    def test_grid_instance_is_initialised( self, mock_Grid_Point, mock_index, mock_volumes ):
        volumes = [ 1.0, 1.0, 1.0, 1.0, 1.0 ]
        mock_volumes.return_value = volumes
        mock_index.side_effect = [ 1, 3 ]
        mock_grid_points = [ Mock( spec=Grid_Point ), Mock( spec=Grid_Point ), Mock( spec=Grid_Point ), Mock( spec=Grid_Point ), Mock( spec=Grid_Point ) ]
        for g in mock_grid_points:
            g.sites = []
        mock_Grid_Point.side_effect = mock_grid_points
        x_coordinates = np.array( [ 0.0, 1.0, 2.0, 3.0, 4.0 ] )
        b = 0.25
        c = 0.1
        limits = [ -0.5, 4.5 ]
        site_set = MagicMock( spec=Set_of_Sites )
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].x = 1.0
        sites[1].x = 3.0
        site_set.__iter__.return_value = iter( sites )
        grid = Grid( x_coordinates=x_coordinates, b=b, c=c, limits=limits, site_set=site_set )
        self.assertEqual( grid.volumes, volumes )
        self.assertEqual( grid.points, mock_grid_points )
        self.x = x_coordinates
        self.limits = limits
        expected_sites_at_grid_points = [ [], [ sites[0] ], [], [ sites[1] ], [] ]
        for p, e in zip( grid.points, expected_sites_at_grid_points ):
            self.assertEqual( p.sites, e )
        # TODO Should really check calls to mocked methods are what we expect
        # TODO The large number of assertions in this test suggests that the Grid __init__ method could be simplifie.

    def test_volumes_from_grid( self ):
        b = 2.0
        c = 3.0
        limits = [-0.5, 4.5]
        grid = np.array( [ 0.0, 1.0, 2.0, 3.0, 4.0 ] )
        with patch( 'project.grid.delta_x_from_grid' ) as mock_delta_x_from_grid:
            mock_delta_x_from_grid.return_value = np.array( [ 1.0, 1.0, 1.0, 1.0, 1.0 ] )
            volumes = volumes_from_grid( b, c, limits, grid )
            expected_volumes = np.array( [ 6.0, 6.0, 6.0, 6.0, 6.0 ] )
            mock_delta_x_from_grid.assert_called_with( grid, limits )
        np.testing.assert_array_equal( volumes, expected_volumes )
       
    def test_delta_x_from_grid( self ):
        grid = np.array( [ 0.0, 1.0, 2.0, 3.0, 4.0 ] )
        limits = [-0.5, 4.5]
        expected_delta_x = np.array( [ 1.0, 1.0, 1.0, 1.0, 1.0 ] )
        delta_x = delta_x_from_grid( grid, limits )
        np.testing.assert_array_equal( expected_delta_x, delta_x )

    def test_delta_x_from_grid_with_uneven_grid( self ):
        grid = np.array( [ 0.5, 1.0, 2.0, 3.0, 5.0 ] )
        limits = [0.0, 5.2]
        expected_delta_x = np.array( [ 0.75, 0.75, 1.0, 1.5, 1.2 ] )
        delta_x = delta_x_from_grid( grid, limits )
        np.testing.assert_array_almost_equal( expected_delta_x, delta_x )
     
if __name__ == '__main__':
    unittest.main()
