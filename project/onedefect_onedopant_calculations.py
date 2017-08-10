from project.defect_species import Defect_Species
from project.site import Site
#from project.practice import *
from project.grid import *
from project.set_of_sites import Set_of_Sites
from project.constants import * 
from project.general_calculations import *
from project.calculation import Calculation

import numpy as np
import scipy
from scipy import stats
import pandas as pd
import math
import itertools
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

from sympy.solvers import solve
from sympy import Symbol, exp, mpmath

from timeit import default_timer as timer

def form_continuum_sites( all_sites, x_min, x_max, n_points, b, c, defect_species ):
    limits = [ x_min, x_max ]    

    grid_1 = np.linspace( x_min, x_max, n_points )
    
    Gd_scaling = len( all_sites.subset( 'Ce' ) ) / len( grid_1 )
    Vo_scaling = len( all_sites.subset( 'O' ) ) / len( grid_1 )
    
    Vo_continuum_grid = Grid( grid_1, b, c, limits, all_sites.subset( 'O' ) )
    Gd_continuum_grid = Grid( grid_1, b, c, limits, all_sites.subset( 'Ce' ) )
    
    Vo_average_energies = np.array( [ site.average_local_energy( method = 'mean' )[0] for site in all_sites.subset( 'O' ) ] )
    Gd_average_energies = np.array( [ site.average_local_energy( method = 'mean' )[0] for site in all_sites.subset( 'Ce' ) ] )
    
    
    Vo_new_energies = griddata( ( [ site.x for site in all_sites.subset( 'O' ) ] ), Vo_average_energies, grid_1, method = 'nearest' )
    Gd_new_energies = griddata( ( [ site.x for site in all_sites.subset( 'Ce' ) ] ), Gd_average_energies, grid_1, method = 'nearest' )
    
    
    Vo_new_sites = Set_of_Sites( [ Site( 'O', x, [ defect_species['Vo'] ], [e], scaling = np.array( Vo_scaling ) ) for x, e in zip( grid_1, Vo_new_energies ) ] )
    Gd_new_sites = Set_of_Sites( [ Site( 'Ce', x, [ defect_species['Gd'] ], [e], scaling = np.array( Gd_scaling ) ) for x, e in zip( grid_1, Gd_new_energies ) ] )   
    
    all_sites_new = Vo_new_sites + Gd_new_sites
    
    return all_sites_new

def example_MF( desired_mobile_defect_mf, slope, intercept ):
    """
    Calculates the mole fraction that is required in the input to get the desired output mole fraction.
    Due to some numerical noise, when the simulation is run the defect mole fractions vary from what they should be and this code uses a solpe and intercept from linear regression of the input and output mole fractions to find what the input mole fraction should be to achieve the desired output.

    Args:
        desired_mobile_defect_MF (float): desired output mole fraction for the mobile defect.
        slope (float): slope calculated from linear regression.
        intercept (float): intercept calcuated from linear regression.

    Returns:
        MF (float): input mole fraction to achieve desired output mole fraction.
    """
    mobile_defect_mf = ( desired_mobile_defect_mf - intercept ) / slope
    MF = [ (mobile_defect_mf), ( mobile_defect_mf ) ]
    return MF


def MF( desired_mobile_defect_mf, slope, intercept ):
    """
    Calculates the mole fraction that is required in the input to get the desired output mole fraction.
    Due to some numerical noise, when the simulation is run the defect mole fractions vary from what they should be and this code uses a solpe and intercept from linear regression of the input and output mole fractions to find what the input mole fraction should be to achieve the desired output.

    Args:
        desired_mobile_defect_MF (float): desired output mole fraction for the mobile defect.
        slope (float): slope calculated from linear regression.
        intercept (float): intercept calcuated from linear regression.

    Returns:
        MF (float): input mole fraction to achieve desired output mole fraction.
    """
    mobile_defect_mf = ( desired_mobile_defect_mf - intercept ) / slope
    MF = [ (mobile_defect_mf), ( mobile_defect_mf * 4 ) ]
    return MF

def calculate_average_molefraction( temp, x_min, x_max, b, c, index, alpha, conv, all_sites, site_labels, boundary_conditions ):
    
    limits = [x_min, x_max ]

    bulk_x_coordinates = np.unique( np.concatenate( ( ( [ x for x in all_sites.get_coords(site_labels[1]) if x <= x_max and x >= x_min ] ), ( [ x for x in all_sites.get_coords(site_labels[0]) if x <= x_max and x >= x_min ] ) ), axis = 0 ) )
    
    bulk_grid = Grid( bulk_x_coordinates, b, c, limits, all_sites )
 
    bulk_mobile_defect_grid = Grid( np.unique( [ x for x in all_sites.get_coords(site_labels[0]) if x <= x_max and x >= x_min ] ), b, c, limits, all_sites.subset( site_labels[0] ) )
    
    phi, rho, niter = calculation( bulk_grid, conv, temp, alpha, boundary_conditions )

    mobile_defect_density = Set_of_Sites( all_sites.subset( site_labels[0] ) ).subgrid_calculate_defect_density( bulk_mobile_defect_grid, bulk_grid, phi, temp )

    mobile_defect_mole_fraction = Set_of_Sites( all_sites.subset( site_labels[0] ) ).calculate_probabilities( bulk_grid, phi, temp )
    mobile_defect_MF = []
    for m in mobile_defect_mole_fraction:
        if m > 0.0:
            mobile_defect_MF.append(m)
    avg_mobile_defect_MF = np.mean(mobile_defect_MF)

    return ( avg_mobile_defect_MF )

#def calculate_GB_properties( temp, x_min, x_max, b, c, index, alpha, conv, desired_mobile_defect_MF, all_sites, site_labels, boundary_conditions ):
    
#    grid = Grid.grid_from_set_of_sites( all_sites, x_min, x_max, b, c)

    
#    calculation_object = Calculation( grid, alpha, conv, temp, boundary_conditions )
#    calculation_object.solve()
#    phi = calculation_object.phi
#    rho = calculation_object.rho
#    niter = calculation_object.niter

#    calculation_object.form_subgrids( site_labels )
#    mobile_defect_grid = calculation_object.subgrids[site_labels[0]]

#    calculation_object.calculate_resistivity_ratio( 0.05 )
#    resistivity_ratio = calculation_object.resistivity_ratio
#    bulk_mobile_defect_density = calculation_object.bulk_mobile_defect_density
   
#    probabilities = all_sites.calculate_probabilities( grid, phi, temp )

#    calculation_object.calculate_mole_fractions()
#    mobile_defect_mole_fraction = calculation_object.mole_fractions[site_labels[0]]
#    dopant_mole_fraction = calculation_object.mole_fractions[site_labels[1] ]

#    return( grid, phi, rho, probabilities, resistivity_ratio, mobile_defect_mole_fraction, dopant_mole_fraction, bulk_mobile_defect_density, niter )
