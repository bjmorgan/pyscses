import numpy as np
from project.variables import wall_potential
import scipy as sp
from scipy.sparse import dia_matrix, diags, spdiags, csc_matrix
from scipy.sparse import linalg 
from scipy.integrate import cumtrapz
from scipy.optimize import brent, minimize_scalar
from project.grid import delta_x_from_grid
from scipy.constants import epsilon_0

def total_charge( c, phi, temp, grid ):
    return sum( grid.charge( phi + c, temp ) )

def total_charge_squared( c, phi, temp, grid ):
    return total_charge( c, phi, temp, grid )**2

class MatrixSolver:
    """ 
    Contains the functions for the finite difference methods used to solve the Poisson-Boltzmann equation on a one-dimensional grid.
    """

    def __init__( self, grid, dielectric, temp, boundary_conditions='dirichlet' ):
        """
        Create a MatrixSolver object.

        Args:
            grid (TODO)
            dielectric ( TODO)
            temp (TODO) Why is this used?
            boundary_conditions (String): Specify the boundary conditions for the matrix solver. 
                                          Allowed values are `dirichlet` and `periodic`.
                                          Default = `dirichlet`.
        Returns:
            None
        """
        self.grid = grid
        self.epsilon = dielectric * epsilon_0
        self.temp = temp 
        allowed_boundary_conditions = [ 'dirichlet', 'periodic' ]
        if boundary_conditions not in allowed_boundary_conditions:
            raise ValueError( boundary_conditions )
        self.boundary_conditions = boundary_conditions
        self.A = self.laplacian_sparse() 

    def laplacian_new_fullmatrix( self ):
        """
        Creates the full tridiagonal matrix used to solve the Poisson-Boltzmann equation using a finite element method.
        deltax is taken as the width of each grid point, from the center point between itself and the next grid point on either side.
        Using a finite difference approximation the diagonal becomes -2.0 / (deltax_1 * deltax_2 ), 
        the upper diagonal becomes 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x2 ) 
        and the lower diagonal becomes 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x1 ).

        Args:
            None

        Returns:
            A (matrix): Full tridiagonal matrix. 
        """
        deltax = self.grid.x[1:] - self.grid.x[:-1] 
        lhs_offset =  (self.grid.x[0] - self.grid.limits[0])
        rhs_offset =  (self.grid.limits[1] - self.grid.x[-1] ) 
        deltax = np.insert( deltax, 0, lhs_offset )
        deltax = np.insert( deltax, len(deltax), rhs_offset )
        delta_x1 = deltax[:-1]
        delta_x2 = deltax[1:]
        diag = -2.0 / (delta_x1 * delta_x2)
        ldiag = 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x1 )
        udiag = 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x2 )
        A = diags( [ diag, udiag[:-1], ldiag[1:] ], [ 0, 1, -1 ] ).A
        if self.boundary_conditions is 'periodic':
            A[0,-1] = ldiag[0]
            A[-1,0] = udiag[-1]
        return A

    def laplacian_sparse( self ):
        """
        Creates a sparse matrix from a full matrix by taking only the tridiagonal values.
        Args:
            None

        Returns:
            A (matrix): Sparse tridiagonal matrix.
        """
        L_sparse = csc_matrix( self.laplacian_new_fullmatrix()) 
        return L_sparse

    def solve( self, phi_old ):
        """ 
        Uses matrix inversion to solve the Poisson-Boltzmann equation through finite difference methods.
    
        Args: 
            phi (array): Electrostatic potential on a one-dimensional grid.
 
        Returns:
            predicted_phi (array): Electrostatic potential on a one-dimensional grid.
        """
        if self.boundary_conditions is 'periodic':
            rho = self.grid.rho( phi_old, self.temp )  
            rho_prime = np.sum( rho * delta_x_from_grid( self.grid.x, self.grid.limits ) ) / np.sum(delta_x_from_grid( self.grid.x, self.grid.limits ) ) 
            rho -= rho_prime
        b = -( rho / self.epsilon )
        predicted_phi = linalg.spsolve( self.A, b )
        return predicted_phi

    

