import numpy as np
import math
from sympy import mpmath
from project.solver import Solver
from project.constants import *
from project.grid import phi_at_x

def sort_xy_data( array ):
    """" Takes an x,y array and sorts it into numerical order from the x axis """
    ind = np.lexsort( ( array[:,1], array[:,0] ) )
    return array[ ind ]

def calculation( grid, conv, temp, alpha ):
    """
    Self-consistent solving of the Poisson-Boltzmann equation. Iterates until the convergence is less than the convergence limit.

    Args:
        grid (object): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
        conv (float): Convergence limit. 
        temp (float): Absolute temperature.
        alpha (float): Damping parameter for the convergence.

    Returns:
        phi (array): Electrostatic potential on a one-dimensional grid. 
        rho (float): Charge density on a one-dimensional grid.
    """
    poisson_solver = Solver( grid.x )

    phi = np.zeros_like( grid.x )
    rho = np.zeros_like( grid.x )
    boundary_condition = poisson_solver.boundary_conditions()

    A = poisson_solver.laplacian_sparse( neumann_bc = [ False, False ] )

    convergence = 1
    while convergence > conv:
        rho = grid.rho( phi, temp )
        predicted_phi = poisson_solver.pb_solver_sparse( rho, boundary_condition, A, dielectric )
        phi =  alpha * predicted_phi + ( 1.0 - alpha ) * phi
        convergence = (sum(( predicted_phi - phi ) **2)) / len( grid.x )
#        print( convergence )
    return phi, rho

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

def r_gb(phi_0, temp, charge ): 
    """
    Mott-Schottky equation for the grain-boundary resistivity ratio
    See e.g. Kim Phys. Chem. Chem. Phys 18, 19787 (2016); Eqn. 8.
    
    Args:
        phi_0  (float): space-charge potential [eV].
        temp   (float): temperature [K].
        charge (float): charge-carrier charge [atomic units]
        
    Returns:
        float: the grain-boundary resistivity ratio [unitless]
    """
    
    return math.exp( charge * phi_0 / ( boltzmann_eV * temp ) ) / ( 2.0 * charge * phi_0 / ( boltzmann_eV * temp ) )

def calculate_activation_energy( ratios, temp ):
    """
    Calculates the ionic conductivity activation energy using the Arrhenius relationship.
    
    Args:
        ratios (array): Grain boundary resistivity at a range of temperatures.
        temp (array): Absolute temperature corresponding to each grain boundary resistivity.

    Returns:
        Ea (float): Ionic conductivity activation energy.
    """
    temp = np.array( temp )
    ratios = np.array( ratios )
    x = 1 / temp
    y = np.log( 1 / ratios )
    slopes = diff_central( x, y )
    Ea = -slopes * boltzmann_eV
    #plt.plot( x, y )
    return Ea

def solve_MS_for_phi(y):
    """
    Solves the Mott-Schottky approximation for the space charge potential usiing the grain bounmdary resistivity.
    
    Args:
        y (float): Grain boundary resistivity

    Returns:
        space charge potential (float): Mott-Schottky approximation.
    """
    return -mpmath.lambertw(-1/(2*y),k=-1)

def debye_length( bulk_density, temp ):
    """
    Calculates the Debye length.
   
    Args: 
        bulk_density (float): Defect density of the mobile defect in the bulk region of the crystal.
        temp (float): Absolute temperature.

    Returns:
        Debye length (float)
    """
 
    return math.sqrt( ( dielectric * vacuum_permittivity * boltzmann_eV * temp ) / ( 2 * ( fundamental_charge ** 2 ) * bulk_density ) )

def space_charge_width( bulk_density, temp, valence, phi_max ):
    """
    Calculates the width of the space charge region.

    Args:
        bulk_density (float): Defect denisty of the mobile defect in the bulk region of the crystal.
        temp (float): Absolute temperature.
        valence (float): Charge of the mobile defect.
        phi_max (float): Space charge potential.

    Returns:
        space charge width (float)
    """
    return 2 * ( debye_length( bulk_density, temp ) * math.sqrt( phi_max / ( ( boltzmann_eV * temp ) / valence ) ) ) 
