import numpy as np
from project.set_of_sites import Set_of_Sites
from project.pbc_solver import PBCSolver
from project.matrix_solver import MatrixSolver
from project.constants import *


class Calculation:

    def __init__(self, grid, alpha, convergence, temp, boundary_conditions):
        self.grid = grid
        self.alpha = alpha
        self.convergence = convergence
        self.temp = temp 
        self.boundary_conditions = boundary_conditions

    def solve( self ):
        """
        Self-consistent solving of the Poisson-Boltzmann equation. Iterates until the convergence is less than the convergence limit.

    Returns:
        phi (array): Electrostatic potential on a one-dimensional grid. 
        rho (float): Charge density on a one-dimensional grid.
        niter (int): Number of iterations performed to reach convergence.
        """

        solvers = { 'dirichlet' : MatrixSolver,
                    'periodic' : PBCSolver }

        poisson_solver = solvers[self.boundary_conditions]( self.grid, dielectric, self.temp )

        phi = np.zeros_like( self.grid.x )
        rho = np.zeros_like( self.grid.x )

        conv = 1
        niter = 0
        while conv > self.convergence:
            rho = self.grid.rho( phi, self.temp )
            predicted_phi = poisson_solver.solve( rho )
            phi =  self.alpha * predicted_phi + ( 1.0 - self.alpha ) * phi
            conv = (sum(( predicted_phi - phi ) **2)) / len( self.grid.x )
            niter += 1

        self.phi = phi
        self.rho = rho
        self.niter = niter

    def form_subgrids( self, site_labels ):
        self.site_labels = site_labels
        subgrids = {} 
        for label in site_labels:
            name = '{}'.format( label ) 
            subgrids[name] = self.grid.subgrid( label )
        self.subgrids = subgrids

    def calculate_resistivity_ratio( self, desired_mobile_defect_MF ):
        mobile_defect_density = Set_of_Sites( self.subgrids[self.site_labels[0]].set_of_sites ).subgrid_calculate_defect_density( self.subgrids[ self.site_labels[0] ], self.grid, self.phi, self.temp )
        
        bulk_mobile_defect_density = len( self.subgrids[self.site_labels[0]].set_of_sites ) * desired_mobile_defect_MF / np.sum( self.subgrids[ self.site_labels[0] ].volumes ) 

        resistivity_ratio = self.subgrids[ self.site_labels[0] ].resistivity_ratio( mobile_defect_density, bulk_mobile_defect_density )
        self.bulk_mobile_defect_density = bulk_mobile_defect_density
        self.resistivity_ratio = resistivity_ratio

    def calculate_mole_fractions( self ):
        mole_fractions = {}
        for label in self.site_labels:
            name = '{}'.format(label)
            mole_fractions[name] = Set_of_Sites(self.subgrids[name].set_of_sites).calculate_probabilities( self.grid, self.phi, self.temp )
        self.mole_fractions = mole_fractions
            
