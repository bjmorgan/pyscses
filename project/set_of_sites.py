import numpy as np
from project.practice import index_of_grid_at_x, phi_at_x
from bisect import bisect_left

def volumes_from_grid( b, c, grid ):
    volumes = np.zeros_like( grid )
    volumes = ( grid[2:] - grid[:-2] ) / 2.0
    volumes = np.insert( volumes, (0, len( volumes ) ), ( grid[1] - grid[0], grid[-1] - grid[-2] ) ) 
    volumes *= ( b * c )
    return volumes
    

class Set_of_Sites:
    def __init__( self, sites ):
        self.sites = sites
    
    def __add__( self, other ):
        if type( other ) is not Set_of_Sites:
            raise TypeError
        return Set_of_Sites( self.sites + other.sites )

    def calculate_rho( self, grid, phi, temp ):
        rho = np.zeros_like( grid )
        charge = np.zeros_like( grid )
        for site in self.sites:
            charge[index_of_grid_at_x( grid, site.x )] += site.charge( phi_at_x( phi, grid, site.x ), temp )
        b = 7.65327e-10
        c = 7.65327e-10 
        # b and c values correct for the 111 2x2 CeO2 grain boundary.
        volumes = volumes_from_grid( b, c, grid )    
        rho = charge / volumes
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

