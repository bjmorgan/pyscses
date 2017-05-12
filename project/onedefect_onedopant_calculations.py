from project.defect_species import Defect_Species
from project.site import Site
from project.solver import Solver
#from project.practice import *
from project.grid import *
from project.set_of_sites import Set_of_Sites
from project.constants import * 
from project.general_calculations import *

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

def site_from_input_file( site, defect_species ):
    """
    Takes the data from the input file and converts it into a site.
    The input data file is a .txt file where each line in the file corresponds to a site.

    Args:
        site (str): A line in the input file.

    Returns:
        Site (object)
    """
    label = site[0]
    x = float(site[1]) 
    defect_labels = site[2::2]
    defect_energies = [ float(e) for e in site[3::2] ]

    return Site( label, x, [ defect_species[l] for l in defect_labels ], defect_energies )


#def form_sites_2( mole_fractions, valence, x_min, x_max, data ):
# def form_sites_2( data_filename, defect_species, x_min, x_max ):

#    GB_data = pd.read_csv( data )
#    Vo_data = GB_data[GB_data.label=='O']
#    Gd_data = GB_data[GB_data.label=='Ce']
#    Vo_label = 'Vo'
#    Gd_label = 'Gd'
#    labels = [ Vo_label, Gd_label ]

#    defect_species = { l : Defect_Species( l, v, m ) for l, v, m in zip( labels, valence, mole_fractions) }

#    Vo_sites = Set_of_Sites( [ Site( 'O', x, [defect_species['Vo']], [e] ) for x, e in zip( Vo_data['x_coord'], Vo_data['seg_energy'] ) if x <= x_max and x >= x_min ] )
#    Gd_sites = Set_of_Sites( [ Site( 'Ce', x, [defect_species['Gd']], [e] ) for x, e in zip( Gd_data['x_coord'], Gd_data['seg_energy'] ) if x <= x_max and x >= x_min ] )
   
#    all_sites = Vo_sites + Gd_sites

 #   new_Vo_sites = all_sites.subset( 'O' )
 #   new_Gd_sites = all_sites.subset( 'Ce' )
 
 #   return( all_sites )
 

def MF( desired_mobile_defect_MF, slope, intercept ):
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
    mobile_defect_MF = ( Desired_mobile_defect_MF - intercept ) / slope
    MF = [ mobile_defect_MF, ( mobile_defect_MF * 4 ) ]
    return MF


def calculate_average_molefraction( temp, x_min, x_max, b, c, index, alpha, conv, all_sites, site_labels ):
    
    bulk_x_coordinates = np.unique( np.concatenate( ( ( [ x for x in all_sites.get_coords(site_labels[1]) if x <= x_max and x >= x_min ] ), ( [ x for x in all_sites.get_coords(site_labels[0]) if x <= x_max and x >= x_min ] ) ), axis = 0 ) )
    
    bulk_grid = Grid( bulk_x_coordinates, b, c, all_sites )
 
    bulk_mobile_defect_grid = Grid( np.unique( [ x for x in all_sites.get_coords(site_labels[0]) if x <= x_max and x >= x_min ] ), b, c, all_sites.subset( site_labels[0] ) )
    
    phi, rho, niter = calculation( bulk_grid, conv, temp, alpha )

    mobile_defect_density = Set_of_Sites( all_sites.subset( site_labels[0] ) ).subgrid_calculate_defect_density( bulk_mobile_defect_grid, bulk_grid, phi, temp )

    mobile_defect_mole_fraction = Set_of_Sites( all_sites.subset( site_labels[0] ) ).calculate_probabilities( bulk_grid, phi, temp )
    mobile_defect_MF = []
    for m in mobile_defect_mole_fraction:
        if m > 0.0:
            mobile_defect_MF.append(m)
    avg_mobile_defect_MF = np.mean(mobile_defect_MF)

    return avg_mobile_defect_MF

def calculate_GB_properties( temp, x_min, x_max, b, c, index, alpha, conv, desired_mobile_defect_MF, all_sites, site_labels ):
    
    x_coordinates = np.unique( np.concatenate( ( ( [ x for x in all_sites.get_coords(site_labels[1]) if x <= x_max and x >= x_min ] ), ( [ x for x in all_sites.get_coords(site_labels[0]) if x <= x_max and x >= x_min ] ) ), axis = 0 ) )

    grid = Grid( x_coordinates, b, c, all_sites )

    mobile_defect_grid = Grid( np.unique( [ x for x in all_sites.get_coords(site_labels[0]) ] ), b, c, all_sites.subset(site_labels[0]) )

    phi, rho, niter = calculation( grid, conv, temp, alpha )

    mobile_defect_density = Set_of_Sites( all_sites.subset(site_labels[0]) ).subgrid_calculate_defect_density( mobile_defect_grid, grid, phi, temp )

    bulk_mobile_defect_density = len( all_sites.subset(site_labels[0]) ) * desired_mobile_defect_MF / np.sum( mobile_defect_grid.volumes )

    resistivity_ratio = mobile_defect_grid.resistivity_ratio( mobile_defect_density, bulk_mobile_defect_density )

    probabilities = all_sites.calculate_probabilities( grid, phi, temp )

    mobile_defect_mole_fraction = Set_of_Sites( all_sites.subset(site_labels[0]) ).calculate_probabilities( grid, phi, temp )
    dopant_mole_fraction = Set_of_Sites( all_sites.subset(site_labels[1]) ).calculate_probabilities( grid, phi, temp )

    return( grid, phi, rho, probabilities, resistivity_ratio, mobile_defect_mole_fraction, dopant_mole_fraction, bulk_mobile_defect_density, niter )

