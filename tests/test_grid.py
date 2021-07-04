import unittest
from pyscses.grid import Grid, delta_x_from_grid
from pyscses.grid_point import GridPoint
from pyscses.set_of_sites import SetOfSites
from pyscses.site import Site
from pyscses.defect_species import DefectSpecies
from unittest.mock import Mock, MagicMock, patch, call
import numpy as np

class TestGrid( unittest.TestCase ):
    @patch('pyscses.grid.index_of_grid_at_x')
    @patch('pyscses.grid.GridPoint')
    def test_grid_instance_is_initialised( self, mock_GridPoint, mock_index ):
        #volumes = [ 0.025, 0.025, 0.025, 0.025, 0.025 ]
        #mock_volumes.return_value = volumes
        mock_index.side_effect = [ 1, 3 ]
        mock_grid_points = [ Mock( spec=GridPoint ), Mock( spec=GridPoint ), Mock( spec=GridPoint ), Mock( spec=GridPoint ), Mock( spec=GridPoint ) ]
        for g in mock_grid_points:
            g.sites = []
        mock_GridPoint.side_effect = mock_grid_points
        x_coordinates = np.array( [ 0.0, 1.0, 2.0, 3.0, 4.0 ] )
        b = 0.25
        c = 0.1
        limits = [ 1.0, 1.0 ]
        limits_for_laplacian = [1.0, 1.0]
        set_of_sites = MagicMock( spec=SetOfSites )
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].x = 1.0
        sites[1].x = 3.0
        sites[0].defect_species = [ Mock( spec=DefectSpecies ) ]
        sites[1].defect_species = [ Mock( spec=DefectSpecies ) ]
        set_of_sites.__iter__.return_value = iter( sites )
        grid = Grid( x_coordinates=x_coordinates, b=b, c=c, limits=limits, limits_for_laplacian=limits_for_laplacian, set_of_sites=set_of_sites )
        self.assertEqual( grid.points, mock_grid_points )
        self.x = x_coordinates
        self.limits = limits
        self.limits_for_laplacian = limits_for_laplacian
        expected_sites_at_grid_points = [ [], [ sites[0] ], [], [ sites[1] ], [] ]
        for p, e in zip( grid.points, expected_sites_at_grid_points ):
            self.assertEqual( p.sites, e )
        for defect_species_list in [ sites[0].defect_species, sites[1].defect_species ]:
            for defect_species in defect_species_list:
                self.assertEqual( defect_species in grid.defect_species, True )
        # TODO Should really check calls to mocked methods are what we expect
        # TODO The large number of assertions in this test suggests that the Grid __init__ method could be simplified / refactored.

    def test_delta_x_from_grid( self ):
        grid = np.array( [ 0.0, 1.0, 2.0, 3.0, 4.0 ] )
        limits = [1.0, 1.0]
        expected_delta_x = np.array( [ 1.0, 1.0, 1.0, 1.0, 1.0 ] )
        delta_x = delta_x_from_grid( grid, limits )
        np.testing.assert_array_equal( expected_delta_x, delta_x )

    def test_delta_x_from_grid_with_uneven_grid( self ):
        grid = np.array( [ 0.5, 1.0, 2.0, 3.0, 5.0 ] )
        limits = [0.5, 1.0]
        expected_delta_x = np.array( [ 0.5, 0.75, 1.0, 1.5, 1.0 ] )
        delta_x = delta_x_from_grid( grid, limits )
        np.testing.assert_array_almost_equal( expected_delta_x, delta_x )

if __name__ == '__main__':
    unittest.main()
