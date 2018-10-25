import numpy as np
from sympy import mpmath
from bisect import bisect_left, bisect_right
from project.set_of_sites import Set_of_Sites
from project.matrix_solver import MatrixSolver
from project.set_up_calculation import calculate_grid_offsets
from project.constants import *
from project.grid import delta_x_from_grid, Grid, phi_at_x
from scipy.optimize import minimize

class Calculation:
    """
    Calculation class
    """

    def __init__(self, grid, bulk_x_min, bulk_x_max, alpha, convergence, dielectric, temp, boundary_conditions, initial_guess=None ):
        self.grid = grid
        self.bulk_x_min = bulk_x_min
        self.bulk_x_max = bulk_x_max
        self.alpha = alpha
        self.convergence = convergence
        self.dielectric = dielectric
        self.temp = temp 
        self.boundary_conditions = boundary_conditions
        self.initial_guess = initial_guess

    def mole_fraction_error( self, input_mole_fractions, target_mole_fractions, approximation):
        """
        Calculates the sum of sqaures error between the output mole fractions calculated from the input mole fractions, compared to the target mole fractions.
        Args:
            input_mole_fractions( list ): Mole fractions for each of the species used in the iterative Poisson-Boltzmann solver.
            target_mole_fractions( list ): The value that the mole fractions should be in the bulk. 

        Returns:
            Sum of squares error between output and target mole fractions. 
        """
        input_mole_fractions = np.array([input_mole_fractions])
        output_mole_fractions = self.mole_fraction_output(input_mole_fractions, approximation)
        squares = []
        for mf in input_mole_fractions:
            for i in range(len(mf)):
                output = output_mole_fractions[0,i]
                target = target_mole_fractions[0,i]
                squares.append(( output - target )**2)
        return sum( squares )

    def mole_fraction_output(self, input_mole_fractions, approximation):
        """
        Calculates the output mole fraction for a given input mole fraction when solving the Poisson-Boltzmann equation.
        Args:
            input_mole_fractions( list ): Mole fractions for each of the species used in the iterative Poisson-Boltzmann solver.
   
        Returns:
            output_mole_fractions( list ): Mole fractions that are calculated from the iterative Poisson-Boltzmann solver.
        """
        for mf in input_mole_fractions:
            for i in range(len(mf)):
                for site in self.grid.set_of_sites.subset( self.site_labels[i] ):
                    for defect in site.defect_species:
                        defect.mole_fraction = input_mole_fractions[0,i]
                    for defect in site.defects:
                        defect.mole_fraction = input_mole_fractions[0,i]
        
        self.solve( approximation )
        species =[]
        for mf in input_mole_fractions:
            for i in range(len(mf)):
                species.append(self.site_labels[i])
        self.form_subgrids( species )
        self.mole_fractions()
        average_mole_fractions = []
        for mf in input_mole_fractions:
            for i in range(len(mf)):
                self.mf[self.site_labels[i]] = [mf for mf in self.mf[self.site_labels[i]] if mf != 0.0 ]
                average_mf = self.calculate_average(self.subgrids[self.site_labels[i]], self.bulk_x_min, self.bulk_x_max, self.mf[self.site_labels[i]])
                average_mole_fractions.append(average_mf) 
        output_mole_fractions = np.array( [ average_mole_fractions ] )
        return output_mole_fractions

        
    def mole_fraction_correction( self, target_mole_fractions, approximation ):
        """
        Starting from an initial guess for the appropriate input mole fractions, minimises the error between the target bulk mole fraction and the output mole fraction from the iterative Poisson-Boltzmann solver. 
        Args:
            target_mole_fractions( list ): The value that the mole fractions should be in the bulk. 

        Returns:
            opt_mole_fractions( list ): The optimum values to be used as the input mole fractions for the iterative Poisson-Boltzmann solver so that the output bulk mole fractions are the target bulk mole fractions. 
        """
        if self.initial_guess == None:
            self.initial_guess = target_mole_fractions
        bounds = [(0.0001, 1)]*len(target_mole_fractions)
        target_mole_fractions = np.array([target_mole_fractions])
        opt_mole_fractions = minimize( self.mole_fraction_error, self.initial_guess, args=(  target_mole_fractions, approximation ), bounds = (bounds) )
        opt_mole_fractions.x = np.array([opt_mole_fractions.x])
        for mf in target_mole_fractions:
            for i in range(len(mf)):
                for site in self.grid.set_of_sites.subset( self.site_labels[i] ):
                    for defect in site.defect_species:
                        defect.mole_fraction = opt_mole_fractions.x[0,i]
                    for defect in site.defects:
                        defect.mole_fraction = opt_mole_fractions.x[0,i]
        self.initial_guess = opt_mole_fractions.x

    def find_index( self, grid, min_cut_off, max_cut_off ):
        min_index = bisect_left( grid.x, min_cut_off )
        max_index = bisect_left( grid.x, max_cut_off )
        return min_index, max_index

    def calculate_offset( self, grid, min_cut_off, max_cut_off ):
        min_index, max_index = self.find_index( grid, min_cut_off, max_cut_off )
        min_offset = ( ( grid.x[ min_index + 1 ] - grid.x[ min_index ] ) / 2 ) + ( ( grid.x[ min_index ] - grid.x[ min_index - 1 ] ) / 2 )
        max_offset = ( ( grid.x[ max_index + 1 ] - grid.x[ max_index ] ) / 2 ) + ( ( grid.x[ max_index ] - grid.x[ max_index - 1 ] ) / 2 )
        return min_offset, max_offset

    def calculate_delta_x( self, grid, min_cut_off, max_cut_off ):
        min_index, max_index = self.find_index( grid, min_cut_off, max_cut_off )
        min_offset, max_offset = self.calculate_offset( grid, min_cut_off, max_cut_off )
        return delta_x_from_grid( grid.x[ min_index+1 : max_index ], [min_offset, max_offset] )

    def calculate_average( self, grid, min_cut_off, max_cut_off, sc_property ):
        min_index, max_index = self.find_index( grid, min_cut_off, max_cut_off )
        delta_x = self.calculate_delta_x( grid, min_cut_off, max_cut_off )
        return np.sum( sc_property[ min_index + 1 : max_index ] * delta_x ) / np.sum( delta_x )

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
        poisson_solver = MatrixSolver( self.grid, self.dielectric, self.temp, boundary_conditions=self.boundary_conditions )

        phi = np.zeros_like( self.grid.x )
        rho = np.zeros_like( self.grid.x )

        conv = 1
        niter = 0
        while conv > self.convergence:
            predicted_phi, rho = poisson_solver.solve( phi )
            if approximation == 'gouy-chapman':
                average_phi = self.calculate_average( self.grid, self.bulk_x_min, self.bulk_x_max, predicted_phi )
                predicted_phi -= average_phi
            if approximation == 'mott-schottky':
                vo_predicted_phi = [ phi_at_x( predicted_phi, self.grid.x, x ) for x in self.grid.subgrid(self.site_labels[0]).x ]
                average_vo_predicted_phi = self.calculate_average( self.grid.subgrid( self.site_labels[0] ), self.bulk_x_min, self.bulk_x_max, predicted_phi )
                predicted_phi -= average_vo_predicted_phi
            phi =  self.alpha * predicted_phi + ( 1.0 - self.alpha ) * phi
            conv = sum( ( predicted_phi - phi )**2) / len( self.grid.x )
            prob = self.grid.set_of_sites.calculate_probabilities( self.grid, phi, self.temp )
            niter += 1
            if niter % 1000 == 0.0:
#            if niter == 1:
                print(conv)
#                print(phi, rho)
#                stop
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
            subgrids[name].volumes[0] = subgrids[name].volumes[1]
            subgrids[name].volumes[-1] = subgrids[name].volumes[1]
        self.subgrids = subgrids

    def create_subregion_sites( self, grid, min_cut_off, max_cut_off ):
        sites = []
        for site in grid.set_of_sites:
            if site.x > min_cut_off and site.x < max_cut_off:
                sites.append( site )
        sites = Set_of_Sites( sites )
        return sites

    def create_space_charge_region( self, grid, pos_or_neg_scr, scr_limit ):
       space_charge_region = []
       self.phi_on_mobile_defect_grid = [ phi_at_x( self.phi, self.grid.x, x ) for x in grid.x ] 
       x_and_phi = np.column_stack( ( grid.x, self.phi_on_mobile_defect_grid ) )
       for i in range( len( x_and_phi ) ):
           if pos_or_neg_scr == 'positive':
               if x_and_phi[i, 1]-x_and_phi[0,1] > scr_limit:
                   space_charge_region.append( x_and_phi[i,0] )
           if pos_or_neg_scr == 'negative': 
               if  x_and_phi[i,1]-x_and_phi[0,1] < scr_limit:
                   space_charge_region.append( x_and_phi[i,0] )
       return space_charge_region

    def calculate_mobile_defect_conductivities( self, pos_or_neg_scr, scr_limit, species, mobility_scaling ):
        space_charge_region = self.create_space_charge_region( self.subgrids[species], pos_or_neg_scr, scr_limit )
        space_charge_region_limits = self.calculate_offset( self.subgrids[species], np.min(space_charge_region), np.max(space_charge_region) )
        space_charge_region_sites = self.create_subregion_sites( self.subgrids[species], np.min(space_charge_region), np.max(space_charge_region) )
        space_charge_region_grid = Grid.grid_from_set_of_sites( space_charge_region_sites, space_charge_region_limits, space_charge_region_limits, self.grid.b, self.grid.c )
        mobile_defect_density = Set_of_Sites( self.subgrids[species].set_of_sites ).subgrid_calculate_defect_density( self.subgrids[species], self.grid, self.phi, self.temp )
        space_charge_region_mobile_defect_density = space_charge_region_sites.subgrid_calculate_defect_density( space_charge_region_grid, self.grid, self.phi, self.temp )
        for site in space_charge_region_sites:
            charge = site.defects[0].valence
            mobilities = site.defects[0].mobility
        if mobility_scaling:
             mobile_defect_conductivity = space_charge_region_mobile_defect_density * ( 1 - space_charge_region_mobile_defect_density ) * charge * mobilities 
        else:
            mobile_defect_conductivity = space_charge_region_mobile_defect_density * charge * mobilities
        min_bulk_index, max_bulk_index = self.find_index( self.subgrids[species], self.bulk_x_min,  self.bulk_x_max )
        self.bulk_limits = self.calculate_offset( self.subgrids[species], self.bulk_x_min, self.bulk_x_max ) 
        bulk_mobile_defect_sites = self.create_subregion_sites( self.subgrids[species], self.bulk_x_min, self.bulk_x_max )
        bulk_mobile_defect_grid = Grid.grid_from_set_of_sites( bulk_mobile_defect_sites, self.bulk_limits, self.bulk_limits, self.grid.b, self.grid.c )
        bulk_mobile_defect_density = Set_of_Sites( bulk_mobile_defect_grid.set_of_sites ).subgrid_calculate_defect_density( bulk_mobile_defect_grid, self.grid, self.phi, self.temp )  
        bulk_mobile_defect_conductivity = bulk_mobile_defect_density * charge * mobilities
        space_charge_array = np.column_stack( ( mobile_defect_conductivity, space_charge_region_grid.x ) )
        bulk_array = np.column_stack( ( bulk_mobile_defect_conductivity, bulk_mobile_defect_grid.x ) )
#        space_charge = sum(mobile_defect_conductivity / space_charge_region_grid.delta_x)
        if mobile_defect_conductivity.all() != 0.0:
            space_charge = sum( space_charge_region_grid.delta_x / mobile_defect_conductivity )
            average_bulk = sum(bulk_mobile_defect_grid.delta_x * bulk_mobile_defect_conductivity) / sum(bulk_mobile_defect_grid.delta_x)
            bulk_1 = [average_bulk] * len( mobile_defect_conductivity)
            bulk = sum(bulk_1 / space_charge_region_grid.delta_x)
#        bulk = sum(bulk_mobile_defect_conductivity / bulk_mobile_defect_grid.delta_x) 
#            ratio = space_charge / bulk        
            ratio = 1 / ( space_charge * bulk )
        else:
            ratio = 0.0 
#        self.depletion_factor = 1 - ( mobile_defect_density / average_bulk )
        return ratio

    def calculate_resistivity_ratio( self, pos_or_neg_scr, scr_limit, mobility_scaling=False ):
        full_conductivity_data = []
        
        for label in self.site_labels:
            c = ( self.calculate_mobile_defect_conductivities( pos_or_neg_scr, scr_limit, label, mobility_scaling  ))
            full_conductivity_data.append(c)
        self.resistivity_ratio = 1 / sum(full_conductivity_data)
            
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
        self.debye_length = math.sqrt( ( self.dielectric * vacuum_permittivity * boltzmann_eV * self.temp ) / ( 2 * ( fundamental_charge ** 2 ) * self.bulk_mobile_defect_density ) )

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
    ratios = np.array(ratios)
    x = 1 / temp
    y = np.log( 1 / ratios )
    slopes = diff_central( x, y )
    Ea = slopes*-boltzmann_eV
    Ea = np.append( Ea, 0 )
    Ea = np.insert( Ea, [0], 0 )
    Ea2 = np.where(Ea == 0, np.nan, Ea )
    return Ea2           
