from copy import copy
from operator import itemgetter
from project.site import Site
import numpy as np
from project.constants import boltzmann_eV
from bisect import bisect_left, bisect_right

def site_from_input_file( site, defect_species, site_charge, core, temperature ):
    """
    Takes the data from the input file and converts it into a site.
    The input data file is a .txt file where each line in the file corresponds to a site.

    Args:
        site (str): A line in the input file.
        defect_species (cls): Class containing information about the defect species present in the system.
        site_charge (bool): True if site charges are to be included in the calculation, false if they are not to be included.

    Returns:
        Site (object)
    """
    label = site[0]
    if site_charge == True:
        valence = float(site[1])
    if site_charge == False:
        valence = 0.0
    x = float(site[2])
    defect_labels = site[3::2]
    defect_energies = [ float(e) for e in site[4::2] ]
    min_energy = min(defect_energies)
    if core == 'single':
        for d_e in defect_energies:
            if d_e > min_energy:
                d_e = 0.0
    if core == 'multi-site':
        for d_e in defect_energies:
            if ( -boltzmann_eV * temperature) <= d_e <= ( boltzmann_eV * temperature ):
                d_e = 0.0
    #defect_energies = [ 0.0 for e in site[4::2] ]
    return Site( label, x, [ defect_species[l] for l in defect_labels ], defect_energies, valence=valence )

def format_line( line, site_charge ):
    """
    Each line in the input file is formatted into separate components to form sites. 
    Args:
        line(str): A line in the input file.
        site_charge (bool): True if site charges are to be included in the calculation, false if they are not to be included.
    Returns: 
        line(list): Formatted line from the input file. 
    """
    # x coordinate
    if site_charge == True:
        line[1] = float( line[1] )
    if site_charge == False:
        line[1] = 0.0

    line[2] = float( line[2] )
    # defect energies
    for i in range( 4, len(line) ):
        line[i] = float( line[i] )
        if line[i] < 0.0:
            line[i] += 0.1
#        line[i] = 0.0
    return line

def load_site_data( filename, x_min, x_max, site_charge ):
    """
    Reads in the input data and formats the input data if the x coordinate values are within the calculation region.
    Args:
        filename(string): Filename for importing the input data.
        x_min(float): Minimum x coordinate value defining the calculation region.
        x_max(float): Maximum x coordinate value defining the calculation region. 
        site_charge (bool): True if site charges are to be included in the calculation, false if they are not to be included.
   
    Return:
        input_data(list): formatted data for calculation.
    """
    with open( filename, 'r' ) as f:
        input_data = [ line.split() for line in f ]
    input_data = [ format_line( line, site_charge ) for line in input_data if ( float(line[2]) > x_min and float(line[2]) < x_max ) ]
    return input_data

def mirror_site_data( site_data, condition = 'symmetrical' ):
    """
    Formatted site data is offset so the maximum x coordinate is shifted to an x coordinate of 0.0. The site data with an x coordinate less than 0.0 is mirrored and the shifted and mirrored site data is concatenated to create a system with two grain boundaries. 
    Args:
        site_data(list): Formatted site data.
    Returns:
        site_data_mirrored(list): Formatted site data for a mirrored system. 
    """

    site_data = sorted(site_data, key=itemgetter(2))
    for line in site_data:
        line[2] = np.round( line[2], 14 )
    if site_data[-1][0] == 'Ce':
        condition = 'symmetrical'
    elif site_data[-1][0] == 'O' and site_data[-5][0] == 'O':
        condition = 'non_symmetrical'
    else:
        raise ValueError('data cannot be mirrored with this end site') 

    if condition == 'symmetrical':
        midpoint = max( [ line[2] for line in site_data ] )
        for line in site_data:
            line[2] -= midpoint
        site_data_mirrored = [ copy(l) for l in site_data if l[2] < 0 ]
        for l in site_data_mirrored:
            l[2] = float(l[2]) * -1 
    if condition == 'non_symmetrical':
        i_list = []
        insert_o_site = copy(site_data[-1])
        offset = site_data[-1][2] - site_data[-5][2]
    
        for i in range( 0, len(site_data) ):
            if round(site_data[i][2], 14) == round(site_data[-1][2], 14):
                i_list.append(i)
        for i in reversed(i_list):
            del site_data[i]
        
        insert_o_site[2] = copy(offset)
        midpoint = max( [ line[2] for line in site_data ] )
        for l in site_data:
            l[2] -= midpoint
    
        site_data_mirrored = [ copy(l) for l in site_data if l[2] < 0 ]
        for l in site_data_mirrored:
            l[2] = float(l[2]) * -1 
        site_data_mirrored = sorted(site_data_mirrored, key=itemgetter(2))
        for l in site_data_mirrored:
            l[2] += offset
        
        site_data_mirrored.insert(0, insert_o_site )
        site_data_mirrored.insert(0, insert_o_site )
        site_data_mirrored.insert(0, insert_o_site )
        site_data_mirrored.insert(0, insert_o_site )

    return site_data + site_data_mirrored

def calculate_grid_offsets( filename, x_min, x_max, system ):
    """
    Reads in the input data calculates the distance to the next site outside of the defined calculation region. Allows calculation of the delta_x and volume values for the endmost grid points.
    Args:
        filename(string): Filename for importing the input data.
        x_min(float): Minimum x coordinate value defining the calculation region.
        x_max(float): Maximum x coordinate value defining the calculation region. 
        system(str): 'single' for single grain boundary systems, 'double' for systems where the input data has been mirrored to give two grain boundaries. 
    Returns:
        limits(list): distance between the midpoint of the endmost sites and the midpoint of the next site outside of the calculation region for the first and last sites respectively. 
        limits_for_laplacian(list): distance between the endmost sites and the next site outside of the calculation region for the first and last sites respectively.
    """
    with open( filename, 'r' ) as f:
        input_data = [ line.split() for line in f ]
        x_coords = np.unique( np.array( [ float(l[2]) for l in input_data ] ) )
        min_index = bisect_left( x_coords, x_min )
        max_index = bisect_right( x_coords, x_max )
        min_offset = ((x_coords[min_index+1]-x_coords[min_index])/2)+((x_coords[min_index]-x_coords[min_index-1])/2)
        max_offset = ((x_coords[max_index+1]-x_coords[max_index])/2)+((x_coords[max_index]-x_coords[max_index-1])/2)
        if system == 'single':
            limits = [ min_offset, max_offset ]
            limits_for_laplacian = [ (x_coords[min_index]-x_coords[min_index-1]), (x_coords[max_index]-x_coords[max_index-1]) ]
        if system == 'double':
            limits = [min_offset, min_offset]
            limits_for_laplacian = [ (x_coords[min_index]-x_coords[min_index-1]), (x_coords[min_index]-x_coords[min_index-1]) ]
    return limits, limits_for_laplacian

