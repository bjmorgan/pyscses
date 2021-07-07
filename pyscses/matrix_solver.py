from __future__ import annotations
import numpy as np
from scipy.sparse import dia_matrix, diags, spdiags, csc_matrix # type: ignore
from scipy.sparse import linalg # type: ignore
from scipy.constants import epsilon_0 # type: ignore
from pyscses.grid import Grid
from typing import Tuple

class MatrixSolver:
    """
    Contains the functions for the finite difference methods used to solve the Poisson-Boltzmann equation on a one-dimensional grid.
    """

    def __init__(self,
                 grid: Grid,
                 dielectric: float,
                 temp: float,
                 boundary_conditions: str = 'dirichlet') -> None:
        """
        Initialise a MatrixSolver object.

        Args:
            grid (Grid): Grid class object which contains the information regarding the calculation grid.
            dielectric (float): dielectric constant for the studied material.
            temp (float): Calculation temperature.
            boundary_conditions (str): Specify the boundary conditions for the matrix solver.
            Allowed values are `dirichlet` and `periodic`. Default = `dirichlet`.

        Returns:
            None

        """
        self.grid = grid
        self.epsilon = dielectric * epsilon_0
        self.temp = temp
        allowed_boundary_conditions = ['dirichlet', 'periodic']
        if boundary_conditions not in allowed_boundary_conditions:
            raise ValueError(boundary_conditions)
        self.boundary_conditions = boundary_conditions
        self.A = self.laplacian_sparse()

    def laplacian_new_fullmatrix(self) -> np.ndarray:
        """
        Creates the full tridiagonal matrix used to solve the Poisson-Boltzmann equation using a finite element method.
        deltax is taken as the width of each grid point, from the center point between itself and the next grid point on either side.
        Using a finite difference approximation the diagonal becomes -2.0 / (deltax_1 * deltax_2 ),
        the upper diagonal becomes 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x2 )
        and the lower diagonal becomes 2.0 / ( ( delta_x1 + delta_x2 ) * delta_x1 ).
        If boundary_conditions are 'periodic', the corner elements of the matrix are filled.

        Args:
            None

        Returns:
            np.ndarray: Full tridiagonal matrix.
        """
        deltax = self.grid.x[1:] - self.grid.x[:-1]
        deltax = np.insert(deltax, len(deltax), self.grid.limits_for_laplacian[0])
        deltax = np.insert(deltax, 0, self.grid.limits_for_laplacian[1])
        delta_x1 = deltax[:-1]
        delta_x2 = deltax[1:]
        diag = -2.0 / (delta_x1 * delta_x2)
        ldiag = 2.0 / ((delta_x1 + delta_x2) * delta_x1)
        udiag = 2.0 / ((delta_x1 + delta_x2) * delta_x2)
        A = diags(diagonals=[diag, udiag[:-1],ldiag[1:]],
                  offsets=[0, 1, -1],
                  format='csc').A
        if self.boundary_conditions == 'periodic':
            A[0,-1] = ldiag[0]
            A[-1,0] = udiag[-1]
        return A

    def laplacian_sparse(self) -> csc_matrix:
        """Creates a sparse matrix from a full matrix by taking only the tridiagonal values.

        Args:
            None

        Returns:
            csc_matrix: Sparse tridiagonal matrix.

        """
        L_sparse = csc_matrix(self.laplacian_new_fullmatrix())
        return L_sparse

    def solve(self,
              phi_old : np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Uses matrix inversion to solve the Poisson-Boltzmann equation through finite difference methods.
        The defect mole fractions are calculated from the electrostatic potential,
        the charge density is calculated from the defect mole fractions and the electrostatic potential
        is then updated using the updated charge density.
        If boundary_conditions are 'periodic', the charge density is minimised before the matrix inversion.
        TODO: Is this really "minimised" for periodic calculations?

        Args:
            phi_old (array): Electrostatic potential on a one-dimensional grid.

        Returns:
            predicted_phi (array): Updated electrostatic potential on a one-dimensional grid.
            rho (array): Updated charge density on a one-dimensional grid.

        """
        rho = self.grid.rho( phi_old, self.temp )
        if self.boundary_conditions == 'periodic':
            rho_prime = np.sum(rho * self.grid.delta_x) / np.sum(self.grid.delta_x)
            rho -= rho_prime
        b = -(rho / self.epsilon)
        predicted_phi = linalg.spsolve(self.A, b)
        return predicted_phi, rho


