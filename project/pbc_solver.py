import numpy as np
from scipy.optimize import minimize_scalar
from project.constants import vacuum_permittivity
from project.grid import delta_x_from_grid


def total_charge( c, phi, temp, grid ):
#    print('total_charge =', sum(grid.charge(phi+c, temp ) ) )
    return sum( grid.charge( phi + c, temp ) )

def total_charge_squared( c, phi, temp, grid ):
    return total_charge( c, phi, temp, grid )**2

def offset( n ):
    return np.insert( n[:-1], 0, n[-1] )

#class PBCSolver:
#    """ Contains the functions for the finite difference methods used to solve the Poisson-Boltzmann equation on a one-dimensional grid. """
#    def __init__( self, grid, dielectric, temp ):
#        self.grid = grid
#        self.epsilon = vacuum_permittivity * dielectric
#        self.temp = temp 
#        self.delta_x = delta_x_from_grid( self.grid.x, self.grid.limits )

#    def integrate( self, n ):
#        n_int = np.cumsum( n * self.delta_x )
#        return n_int - n_int.mean()    

#    def solve( self, rho ):
#        #print('charge density', rho, flush = True)
#        b = -( rho / self.epsilon )
#        predicted_phi = offset( self.integrate( self.integrate(b) ) )
#        #print( 'predicted phi', predicted_phi, flush = True )
#        res = minimize_scalar( total_charge_squared, args=( predicted_phi, self.temp, self.grid ) )
#        corrected_phi = predicted_phi + res.x
#        return corrected_phi

#    def solve( self, phi_old, grid ):
#        res = minimize_scalar( total_charge_squared, args=( phi_old, self.temp, self.grid ) )
#        print(res.x)
#        corrected_phi_old = phi_old + res.x
#        b = -( grid.rho( phi_old + res.x, self.temp ) / self.epsilon )
        #print('b=', b)
#        print('integrate_b =', self.integrate(b) )
#        new_phi = offset( self.integrate( self.integrate(b) ) )
#        return new_phi 


class PBCSolver:
    """
    Poisson solver with periodic boundary conditions.
    Shifts the background potential to minimise the total charge magnitude in the system
    (with rho recalculated), and then performs a double integral to get the
    corresponding potential, phi.
    """

    def __init__( self, grid, dielectric, temp ):
        self.grid = grid
        self.epsilon = vacuum_permittivity * dielectric
        self.temp = temp 
        self.delta_x = delta_x_from_grid( self.grid.x, self.grid.limits )

    def integrate( self, n ):
        n_int = np.roll( np.cumsum( n * self.delta_x ), shift=-1 )
        return n_int - n_int.mean()

#    def solve( self, rho ):
#        b = - rho / self.epsilon
#        phi = self.integrate( self.integrate( b ) )
#        res = minimize_scalar( total_charge_squared, args=( phi, self.temp, self.grid ) )
        # Note: I am still not sure this is entirely consistent. If the potential is shifted to enforce
        # net zero charge, the charge density must also change, and the final potential will 
        # not correspond to the charge density rho fed into this method.
        # A better approach would be to enforce charge neutrality *before* calling PBCSolver.solver, where the
        # shift in phi is added to phi from the *previous* iteration.
#        return phi + res.x


    def shift( self, b ):
        return np.sum( ( np.cumsum( self.delta_x * b ) ) * self.delta_x ) / np.sum( self.delta_x )

    def solve(self, phi_old, grid):
        res = minimize_scalar( total_charge_squared, args=(phi_old, self.temp, self.grid ) )
        corrected_phi_old = phi_old + res.x
        b = -(grid.rho( corrected_phi_old, self.temp ) / self.epsilon )
        new_phi = self.integrate( self.integrate(b) + self.shift(b) )
        return new_phi
        
