from project.site import Site
import numpy as np
from bisect import bisect_left

class Set_of_Sites:
    def __init__( self, grid ):
        self.grid = grid

    def phi_at_x( self, x, phi ):
        index = self.index_of_grid_at_x( x )
        return phi[ index ]

    def index_of_grid_at_x( self, x ):
        return closest_index( self.grid, x )

    def sort_xy_data( array ):
        ind = np.lexsort( ( array[:,1], array[:,0] ) )
        return array[ ind ]

    def calculate_rho( self, sites, phi, temp ):
        rho = np.zeros_like( self.grid )
        for site in sites:
            rho[self.index_of_grid_at_x( site.x )] += site.charge_density( self.phi_at_x( phi, site.x ), temp )
        return rho

    def calculate_defect_density( self, sites, phi, temp ):
        defect_density = np.zeros_like( self.grid )
        for site in sites:
            defect_density[self.index_of_grid_at_x( site.x )] += np.asarray(site.probabilities( self.phi_at_x( phi, site.x ), temp )) * site.density
        return defect_density

def closest_index(myList, myNumber):
    """
    Assumes myList is sorted. Returns index of closest value to myNumber.

    If two numbers are equally close, return the index of the smallest number.
    """
    print( myList, myNumber )
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return len(myList)
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos - 1
