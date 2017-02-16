import numpy as np
from project.variables import wall_potential
from project.constants import vacuum_permittivity
from scipy.sparse import dia_matrix, diags, spdiags, csc_matrix
from scipy.sparse import linalg 

class Solver:
    def __init__( self, grid ):
        self.grid = grid

    def laplacian_new_fullmatrix( self, neumann_bc=[False, False] ):
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
        L_sparse = csc_matrix( self.laplacian_new_fullmatrix( neumann_bc ) )
        return L_sparse

    def boundary_conditions( self ):
        deltax = self.grid[1:] - self.grid[:-1]
        mean_deltax = np.mean( deltax )
        boundary_condition = np.zeros_like( self.grid )
        boundary_condition[ 0 ] = wall_potential * 2.0 / ( mean_deltax * ( mean_deltax + deltax[0] ) )
        boundary_condition[ -1 ] = 0
        #boundary_condition[ -1 ] = other_wall_potential * 2.0 / ( mean_deltax * ( mean_deltax + deltax[-1] ) )
        return boundary_condition

    def epsilon( self, dielectric ):
        return vacuum_permittivity * dielectric

    def pb_solver_fullmatrix( self, rho, boundary_condition, A, dielectric ):
        b = -( rho / self.epsilon( dielectric ) + boundary_condition )
        predicted_phi = np.linalg.solve( A, b )
        return predicted_phi

    def pb_solver_sparse( self, rho, boundary_condition, A, dielectric ):
        b = -( rho / self.epsilon( dielectric ) + boundary_condition )
        predicted_phi = linalg.spsolve( A, b )
        return predicted_phi

