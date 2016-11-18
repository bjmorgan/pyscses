import numpy as np
from project.grid import index_of_grid_at_x, phi_at_x
#from project.grid import volumes_from_grid
from bisect import bisect_left

#def volumes_from_grid( b, c, grid ):
#    volumes = np.zeros_like( grid )
#    volumes = ( grid[2:] - grid[:-2] ) / 2.0
#    volumes = np.insert( volumes, (0, len( volumes ) ), ( grid[1] - grid[0], grid[-1] - grid[-2] ) ) 
#    volumes *= ( b * c )
#    return volumes
    

class Set_of_Sites:
    def __init__( self, sites ):
        self.sites = sites
    
    def __add__( self, other ):
        if type( other ) is not Set_of_Sites:
            raise TypeError
        return Set_of_Sites( self.sites + other.sites )

    def __getitem__( self, index ):
        return self.sites[ index ]

    def calculate_probabilities( self, grid, phi, temp ):
        probability = np.zeros_like( grid.x )
        for site in self.sites:
            probability[index_of_grid_at_x( grid.x, site.x )] = np.asarray( site.probabilities( phi_at_x( phi, grid.x, site.x ), temp ) )
        return probability 

    def calculate_defect_density( self, grid, phi, temp ):
        defect_density = np.zeros_like( grid.x )
        for site in self.sites:
            i = index_of_grid_at_x( grid.x, site.x )
            defect_density[ i ] += np.asarray( site.probabilities( phi_at_x( phi, grid.x, site.x ), temp ) ) / grid.volumes[ i ]
        return defect_density  

