import numpy as np
import math
from project.grid import index_of_grid_at_x, phi_at_x, energy_at_x
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


