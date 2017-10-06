import math
import numpy as np
from sympy import mpmath
import matplotlib.pyplot as plt
from project.set_of_sites import Set_of_Sites
from project.pbc_solver import PBCSolver
from project.matrix_solver import MatrixSolver
from project.general_calculations import diff_central
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
            predicted_phi = poisson_solver.solve( phi, self.grid )
#            predicted_phi = poisson_solver.solve( rho )
            phi =  self.alpha * predicted_phi + ( 1.0 - self.alpha ) * phi
#            plt.plot( self.grid.x, phi )
            conv = (sum(( predicted_phi - phi ) **2)) / len( self.grid.x )
#            print(conv)
#            niter += 1
#            if niter > 0:
#                stop

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

    def solve_MS_approx_for_phi( self, valence ):
        if self.resistivity_ratio < 1.36:
            raise ValueError( "Resistivity ratio < 1.36. Solution not on a real branch." )
        self.ms_phi = (-mpmath.lambertw(-1/(2*self.resistivity_ratio),k=-1)) * ( ( boltzmann_eV * self.temp ) / valence )

    def calculate_debye_length( self ):
        self.debye_length = math.sqrt( ( dielectric * vacuum_permittivity * boltzmann_eV * self.temp ) / ( 2 * ( fundamental_charge ** 2 ) * self.bulk_mobile_defect_density ) )

    def calculate_space_charge_width( self, valence ):
        self.space_charge_width = 2 * ( self.debye_length * math.sqrt( max(self.phi) / ( ( boltzmann_eV * self.temp ) / valence ) ) )


    def mole_fractions( self ):
        mole_fractions = {}
        for label in self.site_labels:
            name = '{}'.format(label)
            mole_fractions[name] = Set_of_Sites(self.subgrids[name].set_of_sites).calculate_probabilities( self.grid, self.phi, self.temp )
        self.mf = mole_fractions
            
