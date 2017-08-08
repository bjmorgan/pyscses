import numpy as np
from scipy.optimize import minimize_scalar
from project.constants import vacuum_permittivity
from project.grid import delta_x_from_grid


def total_charge( c, phi, temp, grid ):
    return sum( grid.rho( phi + c, temp ) )

def total_charge_squared( c, phi, temp, grid ):
    return total_charge( c, phi, temp, grid )**2

def offset( n ):
    return np.insert( n[:-1], 0, n[-1] )

class PBCSolver:
    """ Contains the functions for the finite difference methods used to solve the Poisson-Boltzmann equation on a one-dimensional grid. """
    def __init__( self, grid, dielectric, temp ):
        self.grid = grid
        self.epsilon = vacuum_permittivity * dielectric
        self.temp = temp 

    def integrate( self, n ):
        delta_x = delta_x_from_grid( self.grid.x, self.grid.limits )
        n_int = np.cumsum( n * delta_x )
        return n_int - n_int.mean()    

    def solve( self, rho ):
        b = -( rho / self.epsilon )
        predicted_phi = offset( self.integrate( self.integrate(b) ) )
        res = minimize_scalar( total_charge_squared, args=( predicted_phi, self.temp, self.grid ) )
        corrected_phi = predicted_phi + res.x
        return corrected_phi
 

