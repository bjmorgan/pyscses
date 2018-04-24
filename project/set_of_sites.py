import numpy as np
from project.grid import Grid
from scipy.interpolate import griddata
import math
from project.set_up_calculation import site_from_input_file, load_site_data, mirror_site_data
from project.grid import index_of_grid_at_x, phi_at_x, energy_at_x
from project.site import Site
from project.constants import boltzmann_eV
from bisect import bisect_left

class Set_of_Sites:
    """ The Set_of_Sites object groups together all of the Site objects into one object and contains functions for the calculations that provide properties of all of the sites together rather than individually. """
    def __init__( self, sites ):
        self.sites = sites
    
    def __add__( self, other ):
        """ Allows the concatenation of multiple Set_of_Sites objects"""
        if type( other ) is not Set_of_Sites:
            raise TypeError
        return Set_of_Sites( self.sites + other.sites )

    def __getitem__( self, index ):
        """ Returns the site corresponding to a given index """
        return self.sites[ index ]

    def __iter__( self ):
        """
        Iterator over self.sites
        """
        return iter( self.sites )

    def subset( self, label ):
        """ Returns a list of all the sites which contain a particular defect """
        return [ s for s in self.sites if s.label == label ] 

    def get_coords( self, label ):
        """ Returns a list of the x coordinates for all the sites wich contain a particular defect """
        return [ s.x for s in self.sites if s.label == label ]

    def calculate_energies_on_grid( self, grid, phi ):
        """ 
        Returns an array of energies at their points on a one-dimensional grid.
    
        Args: 
            grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            phi (array): electrostatic potential on a one-dimensional grid. 

        Returns:
            energies_on_grid (array): energies at their grid points
 
        """
        energies_on_grid = np.zeros_like( grid )
        for site in self.sites:
            energies_on_grid[index_of_grid_at_x( grid.x, site.x )] =+ energy_at_x( site.defect_energies, grid.x, site.x )
        return energies_on_grid

    def calculate_probabilities( self, grid, phi, temp ):
        """ 
        Calculates the probability of a site being occupied by its corresponding defect.
    
        Args: 
            grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            phi (array): electrostatic potential on a one-dimensional grid. 
            temp (float): Absolute temperature.

        Returns:
            probability (array): probabilities of defects occupying each site using their grid points
 
        """
        probability = np.zeros_like( grid.x )
        for site in self.sites:
            probability[index_of_grid_at_x( grid.x, site.x )] = np.asarray( site.probabilities( phi_at_x( phi, grid.x, site.x ), temp ) )
        return probability 

    def calculate_defect_density( self, grid, phi, temp ):
        """ 
        Calculates the defect density at each site.
    
        Args: 
            grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            phi (array): electrostatic potential on a one-dimensional grid. 
            temp (float): Absolute temperature.

        Returns:
            defect_density (array): defect density for each site using their grid points
 
        """
        defect_density = np.zeros_like( grid.x )
        for site in self.sites:
            i = index_of_grid_at_x( grid.x, site.x )
            defect_density[ i ] += np.asarray( site.probabilities( phi_at_x( phi, grid.x, site.x ), temp ) ) / grid.volumes[ i ]
        return defect_density  

    def subgrid_calculate_defect_density( self, sub_grid, full_grid, phi, temp ):
        """ 
        Calculates the defect density at each site for a given subset of sites.
    
        Args: 
            subgrid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates. For a given subset of sites. 
            full_grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates. For all sites.
            phi (array): electrostatic potential on a one-dimensional grid. 
            temp (float): Absolute temperature.

        Returns:
            defect_density (array): defect density for each site using their grid points
 
        """
        defect_density = np.zeros_like( sub_grid.x )
        for site in self.sites:
            i = index_of_grid_at_x( sub_grid.x, site.x )
            defect_density[ i ] += np.asarray( site.probabilities( phi_at_x( phi, full_grid.x, site.x ), temp ) ) / sub_grid.volumes[ i ]
        return defect_density  


    def form_continuum_sites( all_sites, x_min, x_max, n_points, b, c, defect_species, limits_for_laplacian ):
        """
        NEEDS UPDATING TO WORK WITH NEW CODE FORMAT.
        INCLUDE DOCSTRING
        """
        limits = [ x_min, x_max ]    

        grid_1 = np.linspace( x_min, x_max, n_points )
    
        Gd_scaling = len( all_sites.subset( 'Ce' ) ) / len( grid_1 )
        Vo_scaling = len( all_sites.subset( 'O' ) ) / len( grid_1 )
    
        Vo_continuum_grid = Grid( grid_1, b, c, limits, limits_for_laplacian, all_sites.subset( 'O' ) )
        Gd_continuum_grid = Grid( grid_1, b, c, limits, limits_for_laplacian, all_sites.subset( 'Ce' ) )
    
        Vo_average_energies = np.array( [ site.average_local_energy( method = 'mean' )[0] for site in all_sites.subset( 'O' ) ] )
        Gd_average_energies = np.array( [ site.average_local_energy( method = 'mean' )[0] for site in all_sites.subset( 'Ce' ) ] )
    
    
        Vo_new_energies = griddata( ( [ site.x for site in all_sites.subset( 'O' ) ] ), Vo_average_energies, grid_1, method = 'nearest' )
        Gd_new_energies = griddata( ( [ site.x for site in all_sites.subset( 'Ce' ) ] ), Gd_average_energies, grid_1, method = 'nearest' )
    
    
        Vo_new_sites = Set_of_Sites( [ Site( 'O', x, [ defect_species['Vo'] ], [e], scaling = np.array( Vo_scaling ) ) for x, e in zip( grid_1, Vo_new_energies ) ] )
        Gd_new_sites = Set_of_Sites( [ Site( 'Ce', x, [ defect_species['Gd'] ], [e], scaling = np.array( Gd_scaling ) ) for x, e in zip( grid_1, Gd_new_energies ) ] )   
    
        all_sites_new = Vo_new_sites + Gd_new_sites
    
        return all_sites_new

    @ classmethod
    def set_of_sites_from_input_data( cls, input_data, limits, defect_species, site_charge ):
        site_data = load_site_data( input_data, limits[0], limits[1], site_charge )
        return Set_of_Sites( [ site_from_input_file( line, defect_species, site_charge ) for line in site_data ] )

    @ classmethod
    def mirrored_set_of_sites_from_input_data( cls, input_data, limits, defect_species ):
        site_data = load_site_data( input_data, limits[0], limits[1] )
        mirrored_data = mirror_site_data( site_data )
        return Set_of_Sites( [ site_from_input_file( line, defect_species ) for line in mirrored_data ] )
