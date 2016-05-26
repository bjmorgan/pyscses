import numpy as np
from project.variables import wall_potential
from project.constants import vacuum_permittivity
from scipy.sparse import dia_matrix, diags
from scipy.sparse import linalg 

class Solver:
    def __init__( self, grid ):
        self.grid = grid

    def laplacian( self ):
        n_points = len( self.grid )
        deltax = self.grid[1:] - self.grid[:-1]
        mean_deltax = np.mean( deltax )
        extended_deltax = np.insert( deltax, (0, len(deltax)), (mean_deltax, mean_deltax) )
        diag = -(extended_deltax[:n_points] + extended_deltax[1:])
        ldiag = extended_deltax[:n_points-1]
        udiag = extended_deltax[2:]
        prefactor = 2.0 / ( extended_deltax[:n_points]* extended_deltax[1:] * (extended_deltax[:n_points] + extended_deltax[1:]) )
        A = diags( [ diag, udiag, ldiag ], [ 0, 1, -1 ] ).A
        L = (A.T * prefactor).T
        return L

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

    def pb_solver( self, rho, boundary_condition, A, dielectric ):
        b = -( rho / self.epsilon( dielectric ) + boundary_condition )
        predicted_phi = np.linalg.solve( A, b )
        return predicted_phi
