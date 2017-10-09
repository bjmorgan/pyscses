import numpy as np
from project.variables import wall_potential
from project.constants import vacuum_permittivity
from scipy.sparse import dia_matrix, diags, spdiags, csc_matrix
from scipy.sparse import linalg 
from scipy.integrate import cumtrapz
from scipy.optimize import minimize_scalar
from project.grid import delta_x_from_grid

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
        self.epsilon = dielectric * vacuum_permittivity
        self.temp = temp # Why is this here?
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
        lhs_offset = self.grid.x[0] - self.grid.limits[0]
        rhs_offset = self.grid.limits[1] - self.grid.x[-1]
        delta_x1 = np.insert( deltax, 0, lhs_offset )
        delta_x2 = np.insert( deltax, len(deltax), rhs_offset )
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
        L_sparse = csc_matrix( self.laplacian_new_fullmatrix( boundary_conditions=boundary_conditions) )
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
            res = minimize_scalar( total_charge_squared, args=(phi_old, self.temp, self.grid ) )
            phi_old += res.x
        b = -( self.grid.rho( phi_old, self.temp ) / self.epsilon )
        predicted_phi = linalg.spsolve( self.A, b )
        return predicted_phi

    

