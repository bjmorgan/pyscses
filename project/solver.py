import numpy as np
from project.variables import wall_potential
from project.constants import vacuum_permittivity
from scipy.sparse import dia_matrix, diags, spdiags, csc_matrix
from scipy.sparse import linalg 

class Solver:
    """ Contains the functions for the finite difference methods used to solve the Poisson-Boltzmann equation on a one-dimensional grid. """
    def __init__( self, grid ):
        self.grid = grid

    def laplacian_new_fullmatrix( self, neumann_bc=[False, False] ):
        """
        Creates the full tridiagonal matrix used to solve the Poisson-Boltzmann equation using a finite element method.
        deltax is taken as the width of each grid point, from the center point between itself and the next grid point on either side.
        Using a finite difference approximation the diagonal becomes -2.0 / (deltax_1 * deltax_2 ), the upper diagonal becomes 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x2 ) and the lower diagonal becomes 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x1 ).

        Args:
            neumann_bc (list) : Default [False, False] for Dirichlet boundary conditions. If True will enforce Neumann boundary conditions.

        Returns:
            A (matrix): Full tridiagonal matrix. 
        """
        n_points = len( self.grid )
        deltax = self.grid[1:] - self.grid[:-1]
        mean_deltax = np.mean( deltax )
        delta_x1 = np.insert( deltax, 0, mean_deltax )
        delta_x2 = np.insert( deltax, len(deltax), mean_deltax )
        diag = -2.0 / (delta_x1 * delta_x2)
        if neumann_bc[0]: # NEEDS TESTING!
            diag[0]  /= 2.0
        if neumann_bc[1]:
            diag[-1] /= 2.0
        ldiag = 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x1 )
        udiag = 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x2 ) 
        A = diags( [ diag, udiag[:-1], ldiag[1:] ], [ 0, 1, -1 ] ).A
        return A

    def laplacian_sparse( self, neumann_bc=[False, False] ):
        """
        Creates a sparse matrix from a full matrix by taking only the tridiagonal values.
        Args:
            neumann_bc (list) : Default [False, False] for Dirichlet boundary conditions. If True will enforce Neumann boundary conditions.

        Returns:
            A (matrix): Sparse tridiagonal matrix.
        """
        L_sparse = csc_matrix( self.laplacian_new_fullmatrix( neumann_bc ) )
        return L_sparse

    def boundary_conditions( self ):
        """ Boundary conditions. Zero for Dirichlet boundary conditions """
        deltax = self.grid[1:] - self.grid[:-1]
        mean_deltax = np.mean( deltax )
        boundary_condition = np.zeros_like( self.grid )
        #boundary_condition[ 0 ] = wall_potential * 2.0 / ( mean_deltax * ( mean_deltax + deltax[0] ) )
        #boundary_condition[ -1 ] = 0
        #boundary_condition[ -1 ] = other_wall_potential * 2.0 / ( mean_deltax * ( mean_deltax + deltax[-1] ) )
        return boundary_condition

    def epsilon( self, dielectric ):
        """ Calculates the vacuum permittivity multiplied by the relative permittivity of a material. """
        return vacuum_permittivity * dielectric

    def pb_solver_fullmatrix( self, rho, boundary_condition, A, dielectric ):
        """
        Uses matrix inversion to solve the Poisson-Boltzmann equation through finite difference methods.
    
        Args: 
            rho (array): Charge denisty on a one-dimensional grid.
            boundary condition (array): Zeros for Dirichlet boundary conditons.
            A (matrix): Full tridiagonal matrix with values from finite difference approximation.
            dielectric (float): Relative permittivity of a material

        Returns:
            predicted_phi (array): Electrostatic potential on a one-dimensional grid.
        """
        b = -( rho / self.epsilon( dielectric ) + boundary_condition )
        predicted_phi = np.linalg.solve( A, b )
        return predicted_phi

    def pb_solver_sparse( self, rho, boundary_condition, A, dielectric ):
        """ 
        Uses matrix inversion to solve the Poisson-Boltzmann equation through finite difference methods.
    
        Args: 
            rho (array): Charge denisty on a one-dimensional grid.
            boundary condition (array): Zeros for Dirichlet boundary conditons.
            A (matrix): Sparse tridiagonal matrix with values from finite difference approximation.
            dielectric (float): Relative permittivity of a material
 
        Returns:
            predicted_phi (array): Electrostatic potential on a one-dimensional grid.
        """

        b = -( rho / self.epsilon( dielectric ) + boundary_condition )
        predicted_phi = linalg.spsolve( A, b )
        return predicted_phi

