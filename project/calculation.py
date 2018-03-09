import math
import numpy as np
from sympy import mpmath
from bisect import bisect_left, bisect_right
from project.set_of_sites import Set_of_Sites
from project.matrix_solver import MatrixSolver
from project.constants import *
from project.grid import delta_x_from_grid, Grid, phi_at_x
from scipy.optimize import minimize

class Calculation:

    def __init__(self, grid, bulk_x_min, bulk_x_max, alpha, convergence, temp, boundary_conditions, initial_guess=None ):
        self.grid = grid
        self.bulk_x_min = bulk_x_min
        self.bulk_x_max = bulk_x_max
        self.alpha = alpha
        self.convergence = convergence
        self.temp = temp 
        self.boundary_conditions = boundary_conditions
        self.initial_guess = initial_guess

    def mole_fraction_error( self, input_mole_fractions, target_mole_fractions):
        """
        Calculates the sum of sqaures error between the output mole fractions calculated from the input mole fractions, compared to the target mole fractions.
        Args:
            input_mole_fractions( list ): Mole fractions for each of the species used in the iterative Poisson-Boltzmann solver.
            target_mole_fractions( list ): The value that the mole fractions should be in the bulk. 

        Returns:
            Sum of squares error between output and target mole fractions. 
        """
        input_mole_fractions = np.array([input_mole_fractions])
        output_mole_fractions = self.mole_fraction_output(input_mole_fractions)
        o_gd = output_mole_fractions[0,1]
        o_vo = output_mole_fractions[0,0]
        t_gd = target_mole_fractions[0,1]
        t_vo = target_mole_fractions[0,0]

        return ( ( o_gd - t_gd )**2 ) + ( ( o_vo - t_vo )**2 )

    def mole_fraction_output(self, input_mole_fractions):
        """
        Calculates the output mole fraction for a given input mole fraction when solving the Poisson-Boltzmann equation.
        Args:
            input_mole_fractions( list ): Mole fractions for each of the species used in the iterative Poisson-Boltzmann solver.
   
        Returns:
            output_mole_fractions( list ): Mole fractions that are calculated from the iterative Poisson-Boltzmann solver.
        """
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
        """
        Starting from an initial guess for the appropriate input mole fractions, minimises the error between the target bulk mole fraction and the output mole fraction from the iterative Poisson-Boltzmann solver. 
        Args:
            target_mole_fractions( list ): The value that the mole fractions should be in the bulk. 

        Returns:
            opt_mole_fractions( list ): The optimum values to be used as the input mole fractions for the iterative Poisson-Boltzmann solver so that the output bulk mole fractions are the target bulk mole fractions. 
        """
        if self.initial_guess == None:
            self.initial_guess = target_mole_fractions
        target_mole_fractions = np.array([target_mole_fractions])
        opt_mole_fractions = minimize( self.mole_fraction_error, self.initial_guess, args=(  target_mole_fractions ), bounds = ( (0.0001,1), (0.0001,1) ) )
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
        self.initial_guess = opt_mole_fractions.x


    def solve( self, approximation ):
        """
        Self-consistent solving of the Poisson-Boltzmann equation. Iterates until the convergence is less than the convergence limit.
        Args:
            approximation( str ): Approximation used for the defect behaviour.
                                  'mott-schottky' - Some defects immobile / fixed to bulk mole fractions.
                                  'gouy-chapman' - All defects mobile / able to redistribute.
        Returns:
            phi (array): Electrostatic potential on a one-dimensional grid. 
            rho (float): Charge density on a one-dimensional grid.
            niter (int): Number of iterations performed to reach convergence.
        """

       
        poisson_solver = MatrixSolver( self.grid, dielectric, self.temp, boundary_conditions=self.boundary_conditions )

        phi = np.zeros_like( self.grid.x )
        rho = np.zeros_like( self.grid.x )

        vo_subgrid = self.grid.subgrid('O')
        min_bulk_index = bisect_left( vo_subgrid.x, self.bulk_x_min )
        max_bulk_index = bisect_left( vo_subgrid.x, self.bulk_x_max )
        conv = 1
        niter = 0
        while conv > self.convergence:
            predicted_phi, rho = poisson_solver.solve( phi )
            if approximation == 'gouy-chapman':
                average_phi = np.sum( predicted_phi[min_bulk_index+1 : max_bulk_index] * delta_x_from_grid( self.grid.x[min_bulk_index+1:max_bulk_index], [self.grid.x[min_bulk_index], self.grid.x[max_bulk_index+1] ] ) ) / np.sum( delta_x_from_grid( self.grid.x[min_bulk_index+1:max_bulk_index], [self.grid.x[min_bulk_index], self.grid.x[max_bulk_index+1] ] ) )  
                predicted_phi -= average_phi
            if approximation == 'mott-schottky':
                vo_predicted_phi = [ phi_at_x( predicted_phi, self.grid.x, x ) for x in vo_subgrid.x ]
                average_vo_predicted_phi = np.sum( vo_predicted_phi[min_bulk_index+1 : max_bulk_index] * delta_x_from_grid( vo_subgrid.x[min_bulk_index+1:max_bulk_index], [vo_subgrid.x[min_bulk_index], vo_subgrid.x[max_bulk_index+1] ] ) ) / np.sum( delta_x_from_grid( vo_subgrid.x[min_bulk_index+1:max_bulk_index], [vo_subgrid.x[min_bulk_index], vo_subgrid.x[max_bulk_index+1] ] ) )  
                predicted_phi -= average_vo_predicted_phi
            phi =  self.alpha * predicted_phi + ( 1.0 - self.alpha ) * phi
            conv = sum( ( predicted_phi - phi )**2) / len( self.grid.x )
            prob = self.grid.set_of_sites.calculate_probabilities( self.grid, phi, self.temp )
            niter += 1
#            if niter % 1000 == 0.0:
            if niter == 1:
                print(phi, rho)
                stop
        self.phi = phi
        self.rho = self.grid.rho( phi, self.temp )
        self.niter = niter

    def form_subgrids( self, site_labels ):
        """
        Creates a Grid class for each species in the system.
        Args:
            site_labels( list ): List of strings for the different site species.
  
        Returns:
            subgrids[name]( cls ): Separate Grid classes for the different site species (name). 
        """
        self.site_labels = site_labels
        subgrids = {} 
        for label in site_labels:
            name = '{}'.format( label ) 
            subgrids[name] = self.grid.subgrid( label )
            subgrids[name].delta_x[0] = subgrids[name].delta_x[1]
            subgrids[name].delta_x[-1] = subgrids[name].delta_x[1]
        self.subgrids = subgrids

    def calculate_resistivity_ratio_new( self ):
        """
        Calculates the ratio of the resistivity in the bulk compared to the resistivity over the space charge regions.
      
        Returns:
             resistivity_ratio( float ): Ratio between resistivity in the bulk and resistivity over the space charge region.  
        """
        scr = []
        mobile_defect_grid = self.subgrids[self.site_labels[0]]
        phi_at_mdg = [ phi_at_x( self.phi, self.grid.x, x ) for x in mobile_defect_grid.x ]
        x_and_phi = np.column_stack(( mobile_defect_grid.x, phi_at_mdg ))
        for i in range( len( x_and_phi ) ):
            if x_and_phi[i, 1]-x_and_phi[0,1] > 2e-2:
                scr.append( x_and_phi[i,0] )
        x_min_scr = min( scr )
        x_max_scr = max( scr )
        min_scr_index = bisect_left( mobile_defect_grid.x, x_min_scr )
        max_scr_index = bisect_left( mobile_defect_grid.x, x_max_scr )
        min_scr_offset = ((mobile_defect_grid.x[min_scr_index+2]-mobile_defect_grid.x[min_scr_index+1])/2) + ((mobile_defect_grid.x[min_scr_index+1]-mobile_defect_grid.x[min_scr_index])/2)
        max_scr_offset = ((mobile_defect_grid.x[max_scr_index+1]-mobile_defect_grid.x[max_scr_index])/2) + ((mobile_defect_grid.x[max_scr_index]-mobile_defect_grid.x[max_scr_index-1])/2)
        scr_limits = [ min_scr_offset, max_scr_offset ]
        scr_sites = []
        for site in mobile_defect_grid.set_of_sites:
            if site.x > x_min_scr and site.x < x_max_scr:
                scr_sites.append(site)
        scr_set_of_sites = Set_of_Sites(scr_sites)
        scr_grid = Grid.grid_from_set_of_sites( scr_set_of_sites, scr_limits, scr_limits, self.grid.b, self.grid.c )
        mobile_defect_density = scr_set_of_sites.subgrid_calculate_defect_density( scr_grid, self.grid, self.phi, self.temp )
        min_bulk_index = bisect_left( mobile_defect_grid.x, self.bulk_x_min )
        max_bulk_index = bisect_left( mobile_defect_grid.x, self.bulk_x_max )
        min_bulk_offset = ((mobile_defect_grid.x[min_bulk_index+2]-mobile_defect_grid.x[min_bulk_index+1])/2) + ((mobile_defect_grid.x[min_bulk_index+1]-mobile_defect_grid.x[min_bulk_index])/2)
        max_bulk_offset = ((mobile_defect_grid.x[max_bulk_index+1]-mobile_defect_grid.x[max_bulk_index])/2) + ((mobile_defect_grid.x[max_bulk_index]-mobile_defect_grid.x[max_bulk_index-1])/2)
        self.bulk_limits = [ min_bulk_offset, max_bulk_offset ]
        bulk_mobile_defect_sites = []
        for site in mobile_defect_grid.set_of_sites:
            if site.x > self.bulk_x_min and site.x < self.bulk_x_max:
                bulk_mobile_defect_sites.append(site)
        bulk_mobile_defect_set_of_sites = Set_of_Sites(bulk_mobile_defect_sites)
        bulk_mobile_defect_grid = Grid.grid_from_set_of_sites( bulk_mobile_defect_set_of_sites, self.bulk_limits, self.bulk_limits, self.grid.b, self.grid.c )
        bulk_mobile_defect_density = Set_of_Sites( bulk_mobile_defect_grid.set_of_sites ).subgrid_calculate_defect_density( bulk_mobile_defect_grid, self.grid, self.phi, self.temp )
        self.resistivity_ratio = sum( bulk_mobile_defect_density / bulk_mobile_defect_grid.delta_x ) * sum( scr_grid.delta_x / mobile_defect_density )        
    

    def solve_MS_approx_for_phi( self, valence ):
        """
        Using the calculated resistivity ratio, solves the Mott-Schottky model to calculate the space charge potential (analogous to experimental analysis).
        If the resistivity ratio calculated is less that 1.36, ValueError raised as solution to the Mott-Schottky model is not on a real branch.
        Args:
            valence( float ): Charge of species at site. 
   
        Returns:
            ms_phi( float ): Space charge potential calculated from Mott-Schottky model
        """
        if self.resistivity_ratio < 1.36:
            raise ValueError( "Resistivity ratio < 1.36. Solution not on a real branch." )
        self.ms_phi = (-mpmath.lambertw(-1/(2*self.resistivity_ratio),k=-1)) * ( ( boltzmann_eV * self.temp ) / valence )

    def calculate_debye_length( self ):
        """
        Returns:
            debye_length( float ): Debye length as derived from Poisson-Boltzmann equation
        """
        self.debye_length = math.sqrt( ( dielectric * vacuum_permittivity * boltzmann_eV * self.temp ) / ( 2 * ( fundamental_charge ** 2 ) * self.bulk_mobile_defect_density ) )

    def calculate_space_charge_width( self, valence ):
        """
        Calculates the approximate space charge width from the debye length.
   
        Returns:
            space_charge_width( float ): Approximate space charge width.
        """
        self.space_charge_width = 2 * ( self.debye_length * math.sqrt( max(self.phi) / ( ( boltzmann_eV * self.temp ) / valence ) ) )


    def mole_fractions( self ):
        """
        Calculates the mole fractions (probability of defects occupation ) for each site on the subgrid for each species.
 
        Returns:
            mf( np.array ): defect species mole fractions for each site on the subgrid for each site species. 
        """
        mole_fractions = {}
        for label in self.site_labels:
            name = '{}'.format(label)
            mole_fractions[name] = Set_of_Sites(self.subgrids[name].set_of_sites).calculate_probabilities( self.grid, self.phi, self.temp )
        self.mf = mole_fractions


def diff_central(x, y):
    """
    Calculates the numerical derivative of x,y data using a central difference approximation.
   
    Args:
        x, y (array): x,y data to be differentiated.     
    Returns:
        (array): numerical derivative of x,y
    """
    x0 = x[:-2]
    x1 = x[1:-1]
    x2 = x[2:]
    y0 = y[:-2]
    y1 = y[1:-1]
    y2 = y[2:]
    f = (x2 - x1)/(x2 - x0)
    return (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0) 

def calculate_activation_energies( ratios, temp ):
    """
    Solves the Arrhenius equation using the calculated resistivity ratio for a series of temperatures to calculate the activation energy for single defect segregation.
    Uses a central difference approach so endpoints filled in with NaN to create an array with the same length as the args.
   
    Args:
        ratios( list ): Resistivity ratios calculated for a range of temperatures.
        temp( list ): Temperature (values used to calculate resistivity ratio values).

    Returns:
        Ea2( np.array ): The activation energy calculated for different temperatures.  
    """
    temp = np.array( temp )
    x = 1 / temp
    y = np.log( 1 / ratios )
    slopes = diff_central( x, y )
    Ea = slopes*-boltzmann_eV
    Ea = np.append( Ea, 0 )
    Ea = np.insert( Ea, [0], 0 )
    Ea2 = np.where(Ea == 0, np.nan, Ea )
    return Ea2           
