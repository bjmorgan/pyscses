from project.defect_species import Defect_Species
from project.defect_at_site import Defect_at_Site
from project.site import Site
from project.solver import Solver
from project.grid import *
from project.set_of_sites import Set_of_Sites
from project.constants import boltzmann_eV, Gd_doped_ceria_dielectric as dielectric

import numpy as np
import math
import itertools
import matplotlib
import matplotlib.pyplot as plt

from timeit import default_timer as timer
from scipy.interpolate import griddata

def calculation( grid, conv, temp, alpha ):
    poisson_solver = Solver( grid.x )
    
    phi = np.zeros_like( grid.x )
    rho = np.zeros_like( grid.x )
    boundary_condition = poisson_solver.boundary_conditions()
    
    A = poisson_solver.laplacian_sparse( neumann_bc = [ False, True ] )
    
    convergence = 1
    while convergence > conv:
        rho = grid.rho( phi, temp ) 
        predicted_phi = poisson_solver.pb_solver_sparse( rho, boundary_condition, A, dielectric )
        phi =  alpha * predicted_phi + ( 1.0 - alpha ) * phi
        convergence = (sum(( predicted_phi - phi ) **2)) / len( grid.x )
        #print(convergence)
    return phi, rho

def output( grid, phi, rho, temp, Vo_sites, Gd_sites, index, m, n, conv, fixed, Gd_dist_temp, data ):
    defect_density_plus = Vo_sites.calculate_defect_density( grid, phi, temp )
    defect_density_minus = Gd_sites.calculate_defect_density( grid, phi, temp )  
    site_set = Vo_sites + Gd_sites
        
    probabilities = site_set.calculate_probabilities( grid, phi, temp )      
    plt.plot( grid.x, phi, color = "blue", label = 'Potential' ) 
    plt.show()
    plt.plot( grid.x, rho )
    plt.show()
    Vo_mole_fraction = Vo_sites.calculate_probabilities( grid, phi, temp )
    Gd_mole_fraction = Gd_sites.calculate_probabilities( grid, phi, temp )

    plt.plot( grid.x, Vo_mole_fraction, color = 'red' )
    plt.show()
    plt.plot( grid.x, Gd_mole_fraction, color = 'pink' )
    plt.show()
    print( 'max phi =', max( phi ) )

    file_name_1 = '{}_continuum_data_molefracs{}_conv{}_npoints{}_temp{}_fixed{}_fixedtemp{}.dat'.format( index, m, conv, n, temp, fixed, Gd_dist_temp )
    header = '{}_continuum_data_molefracs{}_conv{}_temp{}_fixed{}_fixedtemp{}_maxpotential{}'"\n x phi rho".format( index, m, conv, temp, fixed, Gd_dist_temp, max(phi) )    
    np.savetxt( file_name_1, data, header = header )

