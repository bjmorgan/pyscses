import numpy as np
from project.practice import index_of_grid_at_x, phi_at_x
from bisect import bisect_left

class Set_of_Sites:
    def __init__( self, sites ):
        self.sites = sites
    
    def __add__( self, other ):
        if type( other ) is not Set_of_Sites:
            raise TypeError
        return Set_of_Sites( self.sites + other.sites )

    def calculate_rho( self, grid, phi, temp ):
        rho = np.zeros_like( grid )
        for site in self.sites:
            rho[index_of_grid_at_x( grid, site.x )] += site.charge_density( phi_at_x( phi, grid, site.x ), temp )
        return rho 

    def calculate_probabilities( self, grid, phi, temp ):
        probability = np.zeros_like( grid )
        for site in self.sites:
            probability[index_of_grid_at_x( grid, site.x )] = np.asarray( site.probabilities( phi_at_x( phi, grid, site.x ), temp ) )
        return probability 

    def calculate_defect_density( self, grid, phi, temp ):
        defect_density = np.zeros_like( grid )
        for site in self.sites:
            defect_density[index_of_grid_at_x( grid, site.x )] += np.asarray(site.probabilities( phi_at_x( phi, grid, site.x ), temp )) * site.density
        return defect_density  

