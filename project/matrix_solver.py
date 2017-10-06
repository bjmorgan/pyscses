import numpy as np
from project.variables import wall_potential
from project.constants import vacuum_permittivity
from scipy.sparse import dia_matrix, diags, spdiags, csc_matrix
from scipy.sparse import linalg 
from scipy.integrate import cumtrapz
from scipy.optimize import minimize_scalar
from project.grid import delta_x_from_grid
from project.pbc_solver import total_charge_squared

class MatrixSolver:
    """ Contains the functions for the finite difference methods used to solve the Poisson-Boltzmann equation on a one-dimensional grid. """
    def __init__( self, grid, dielectric, temp ):
        self.grid = grid
        self.epsilon = dielectric * vacuum_permittivity
        self.temp = temp
        self.A = self.laplacian_sparse() 

    def laplacian_new_fullmatrix( self ):
        """
        Creates the full tridiagonal matrix used to solve the Poisson-Boltzmann equation using a finite element method.
        deltax is taken as the width of each grid point, from the center point between itself and the next grid point on either side.
        Using a finite difference approximation the diagonal becomes -2.0 / (deltax_1 * deltax_2 ), the upper diagonal becomes 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x2 ) and the lower diagonal becomes 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x1 ).

        Args:
            None

        Returns:
            A (matrix): Full tridiagonal matrix. 
        """
        deltax = self.grid.x[1:] - self.grid.x[:-1]
        mean_deltax = np.mean( deltax )
        delta_x1 = np.insert( deltax, 0, mean_deltax )
        delta_x2 = np.insert( deltax, len(deltax), mean_deltax )
        diag = -2.0 / (delta_x1 * delta_x2)
        ldiag = 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x1 )
        udiag = 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x2 ) 
        A = diags( [ diag, udiag[:-1], ldiag[1:] ], [ 0, 1, -1 ] ).A
        A[0,-1] = A[0,1]
        A[-1,0] = A[-1, -2]
        return A

    def laplacian_sparse( self ):
        """
        Creates a sparse matrix from a full matrix by taking only the tridiagonal values.
        Args:
            None

        Returns:
            A (matrix): Sparse tridiagonal matrix.
        """
        L_sparse = csc_matrix( self.laplacian_new_fullmatrix() )
        return L_sparse

#    def solve( self, rho ):
#        """ 
#        Uses matrix inversion to solve the Poisson-Boltzmann equation through finite difference methods.
#    
#        Args: 
#            rho (array): Charge denisty on a one-dimensional grid.
# 
#        Returns:
#            predicted_phi (array): Electrostatic potential on a one-dimensional grid.
#        """
#        b = -( rho / self.epsilon )
#        predicted_phi = linalg.spsolve( self.A, b )
#        return predicted_phi


    def integrate( self, n ):
        n_int = np.roll( np.cumsum( n * self.delta_x ), shift=-1 )
        return n_int - n_int.mean()

    def shift( self, b ):
        return np.sum( ( np.cumsum( self.delta_x * b ) ) * self.delta_x ) / np.sum( self.delta_x )

    def solve(self, phi_old, grid):
        res = minimize_scalar( total_charge_squared, args=(phi_old, self.temp, self.grid ) )
        corrected_phi_old = phi_old + res.x
        b = -(grid.rho( corrected_phi_old, self.temp ) / self.epsilon )
        new_phi = linalg.spsolve( self.A, b )
        return new_phi

    

#    def laplacian_new_fullmatrix( self ):
#        deltax = self.grid.x[1:] - self.grid.x[:-1]
#        mean_deltax = np.mean( deltax )
#        delta_x1 = np.insert( deltax, 0, mean_deltax )
#        delta_x2 = np.insert( deltax, len(deltax), mean_deltax )
#        diag = -2.0 / (delta_x1 * delta_x2)
#        ldiag = 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x1 )
#        udiag = 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x2 )
#        #A = diags( [ diag, udiag[:-1], ldiag[1:] ], [ 0, 1, -1 ] ).A
#        A = np.diag( diag, 0 ) + np.diag( ldiag, -1 ) + np.diag( udiag, +1 )
#        A[0,len(diag-1)] = A[0,len(diag-1)]+1
#        A[len(diag-1),0] = A[len(diag-1),0]+1
#        return A

#    def solve(self, rho):
#        b = -( rho / self.epsilon )
#        predicted_phi = np.zeros_like(rho)
#        predicted_phi[1:] = linalg.spsolve(A[1:,1:],b[1:]) 
#        predicted_phi[0] = predicted_phi[-1]
#        return predicted_phi
