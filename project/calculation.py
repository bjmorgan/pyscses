import math
import numpy as np
from sympy import mpmath
import matplotlib.pyplot as plt
from bisect import bisect_left, bisect_right
from project.set_of_sites import Set_of_Sites
from project.pbc_solver import PBCSolver
from project.matrix_solver import MatrixSolver
from project.general_calculations import diff_central
from project.constants import *
from project.grid import delta_x_from_grid, Grid
from scipy.optimize import minimize

class Calculation:

    def __init__(self, grid, bulk_x_min, bulk_x_max, alpha, convergence, temp, boundary_conditions ):
        self.grid = grid
        self.bulk_x_min = bulk_x_min
        self.bulk_x_max = bulk_x_max
        self.alpha = alpha
        self.convergence = convergence
        self.temp = temp 
        self.boundary_conditions = boundary_conditions

    def mole_fraction_error( self, input_mole_fractions, target_mole_fractions):
        input_mole_fractions = np.array([input_mole_fractions])
        output_mole_fractions = self.mole_fraction_output(input_mole_fractions)
        o_gd = output_mole_fractions[0,1]
        o_vo = output_mole_fractions[0,0]
        t_gd = target_mole_fractions[0,1]
        t_vo = target_mole_fractions[0,0]

        return ( ( o_gd - t_gd )**2 ) + ( ( o_vo - t_vo )**2 )

    def mole_fraction_output(self, input_mole_fractions):
        for site in self.grid.set_of_sites.subset('O'):
            for defect in site.defect_species:
                defect.mole_fraction = input_mole_fractions[0,0]
            for defect in site.defects:
                defect.mole_fraction = input_mole_fractions[0,0]
        for site in self.grid.set_of_sites.subset('Ce'):
            for defect in site.defect_species:
                defect.mole_fraction = input_mole_fractions[0,1]
            for defect in site.defects:
                defect.mole_fraction = input_mole_fractions[0,1]
        Vo_molfracs = []
        Gd_molfracs = []
        self.solve()
        self.form_subgrids( [ 'O', 'Ce' ] )
        self.mole_fractions()
        min_vo_index = bisect_left( self.subgrids['O'].x, self.bulk_x_min )
        max_vo_index = bisect_left( self.subgrids['O'].x, self.bulk_x_max )
        min_gd_index = bisect_left( self.subgrids['Ce'].x, self.bulk_x_min )
        max_gd_index = bisect_left( self.subgrids['Ce'].x, self.bulk_x_max )
        self.mf['O'] = [ mf for mf in self.mf['O'] if mf != 0.0 ]
        self.mf['Ce'] = [ mf for mf in self.mf['Ce'] if mf != 0.0 ]

        for mf in self.mf['O'][min_vo_index+1 : max_vo_index]:
            if mf > 0.0:
                Vo_molfracs.append(mf)

        avg_vo_molfrac = np.sum( Vo_molfracs * delta_x_from_grid( self.subgrids['O'].x[min_vo_index+1:max_vo_index], [ self.subgrids['O'].x[min_vo_index], self.subgrids['O'].x[max_vo_index+1] ] ) ) / np.sum( delta_x_from_grid( self.subgrids['O'].x[min_vo_index+1:max_vo_index], [ self.subgrids['O'].x[min_vo_index], self.subgrids['O'].x[max_vo_index+1] ] ) ) 
       
        for mf in self.mf['Ce'][min_gd_index+1 : max_gd_index]:
            if mf > 0.0:
                Gd_molfracs.append(mf)

        avg_gd_molfrac = np.sum( Gd_molfracs * delta_x_from_grid( self.subgrids['Ce'].x[min_gd_index+1:max_gd_index], [ self.subgrids['Ce'].x[min_gd_index], self.subgrids['Ce'].x[max_gd_index+1] ] ) ) / np.sum( delta_x_from_grid( self.subgrids['Ce'].x[min_gd_index+1:max_gd_index], [ self.subgrids['Ce'].x[min_gd_index], self.subgrids['Ce'].x[max_gd_index+1] ] ) ) 
        
        output_mole_fractions = np.array( [ [ avg_vo_molfrac, avg_gd_molfrac ] ] )
        return output_mole_fractions

        
    def mole_fraction_correction( self, target_mole_fractions ):
        target_mole_fractions = np.array([target_mole_fractions])
        initial_guess = np.array( [ 0.05, 0.2 ] )
        opt_mole_fractions = minimize( self.mole_fraction_error, initial_guess, args=(  target_mole_fractions ), bounds = ( (0.0001,1), (0.0001,1) ) )
        opt_mole_fractions.x = np.array([opt_mole_fractions.x])
        for site in self.grid.set_of_sites.subset('O'):
            for defect in site.defect_species:
                defect.mole_fraction = opt_mole_fractions.x[0,0]
            for defect in site.defects:
                defect.mole_fraction = opt_mole_fractions.x[0,0]
        for site in self.grid.set_of_sites.subset('Ce'):
            for defect in site.defect_species:
                defect.mole_fraction = opt_mole_fractions.x[0,1]
            for defect in site.defects:
                defect.mole_fraction = opt_mole_fractions.x[0,1]


    def form_bulk_grid( self ):
        self.min_bulk_index = bisect_left( self.grid.x, self.bulk_x_min )
        self.max_bulk_index = bisect_left( self.grid.x, self.bulk_x_max )
        self.bulk_limits = [ self.grid.x[self.min_bulk_index], self.grid.x[self.max_bulk_index+1 ] ]
        bulk_sites = []
        for site in self.grid.set_of_sites:
            if site.x > self.bulk_x_min and site.x < self.bulk_x_max:
                bulk_sites.append(site)
        bulk_set_of_sites = Set_of_Sites(bulk_sites)
        bulk_grid = Grid.grid_from_set_of_sites( bulk_set_of_sites, self.bulk_limits[0], self.bulk_limits[1], self.grid.b, self.grid.c )    
        return bulk_grid

    def solve( self ):
        """
        Self-consistent solving of the Poisson-Boltzmann equation. Iterates until the convergence is less than the convergence limit.

    Returns:
        phi (array): Electrostatic potential on a one-dimensional grid. 
        rho (float): Charge density on a one-dimensional grid.
        niter (int): Number of iterations performed to reach convergence.
        """

        poisson_solver = MatrixSolver( self.grid, dielectric, self.temp, boundary_conditions=self.boundary_conditions )

        phi = np.zeros_like( self.grid.x )
        rho = np.zeros_like( self.grid.x )

        min_bulk_index = bisect_left( self.grid.x, self.bulk_x_min )
        max_bulk_index = bisect_left( self.grid.x, self.bulk_x_max )
        conv = 1
        niter = 0
        while conv > self.convergence:
            predicted_phi = poisson_solver.solve( phi )
            average_phi = np.sum( predicted_phi[min_bulk_index+1 : max_bulk_index] * delta_x_from_grid( self.grid.x[min_bulk_index+1:max_bulk_index], [self.grid.x[min_bulk_index], self.grid.x[max_bulk_index+1] ] ) ) / np.sum( delta_x_from_grid( self.grid.x[min_bulk_index+1:max_bulk_index], [self.grid.x[min_bulk_index], self.grid.x[max_bulk_index+1] ] ) )  
#            predicted_phi -=  predicted_phi[0]
            predicted_phi -= average_phi
            phi =  self.alpha * predicted_phi + ( 1.0 - self.alpha ) * phi
            conv = sum( ( predicted_phi - phi )**2) / len( self.grid.x )
            niter += 1
            if niter % 5000== 0:
                print(conv, flush=True)
        self.phi = phi
        self.average_phi = np.sum( predicted_phi[min_bulk_index+1 : max_bulk_index] * delta_x_from_grid( self.grid.x[min_bulk_index+1:max_bulk_index], [self.grid.x[min_bulk_index], self.grid.x[max_bulk_index+1] ] ) ) / np.sum( delta_x_from_grid( self.grid.x[min_bulk_index+1:max_bulk_index], [self.grid.x[min_bulk_index], self.grid.x[max_bulk_index+1] ] ) )  
        self.rho = self.grid.rho( phi, self.temp )
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
            
