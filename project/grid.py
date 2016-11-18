import numpy as np
from bisect import bisect_left

def phi_at_x( phi, grid, x ):
    index = index_of_grid_at_x( grid, x )
    return phi[ index ]

def index_of_grid_at_x( grid, x ):
    return closest_index( grid, x )
#     return bisect( grid, x )

def closest_index(myList, myNumber):
    """
    Assumes myList is sorted. Returns index of closest value to myNumber.

    If two numbers are equally close, return the index of the smallest number.
    """
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

def volumes_from_grid( b, c, grid ):
    volumes = np.zeros_like( grid )
    volumes = ( grid[2:] - grid[:-2] ) / 2.0
    volumes = np.insert( volumes, (0, len( volumes ) ), ( grid[1] - grid[0], grid[-1] - grid[-2] ) )
    volumes *= ( b * c )
    return volumes

class Grid_Point:

  def __init__( self, x, volume ):
      self.x = x
      self.volume = volume
      self.sites = []

class Grid:

  def __init__( self, x_coordinates, b, c, site_set ):
      # x_coordinates need to be sorted for the delta_x calculation in volumes_from_grid
      self.volumes = volumes_from_grid( b, c, x_coordinates )
      self.points = [ Grid_Point( x, v ) for x, v in zip( x_coordinates, self.volumes ) ]
      self.x = x_coordinates
      for site in site_set:
          i = index_of_grid_at_x( self.x, site.x )
          self.points[ i ].sites.append( site ) 

  def __getitem__( self, key ):
      return self.points[ key ]

  def charge( self, phi, temp ):
      rho = np.zeros_like( self.x )
      charge = np.zeros_like( self.x )
      for i, point in enumerate( self ):
          charge[ i ] = sum( [ site.charge( phi_at_x( phi, self.x, site.x ), temp ) for site in point.sites ] )
      return charge

  def rho( self, phi, temp ):
      return self.charge( phi, temp ) / self.volumes


#x_coordinates = np.array( [1.0, 2.0, 3.0, 4.0] )
#b = 7.65327e-10
#c = 7.65327e-10
#g = Grid( x_coordinates, b, c )

#print( g[1].x )
#print( g.x )

