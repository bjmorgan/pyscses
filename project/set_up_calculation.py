from project.site import Site
import numpy as np
from bisect import bisect_left, bisect_right

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
    valence = float(site[1])
#    valence = 0.0
    x = float(site[2])
    defect_labels = site[3::2]
    defect_energies = [ float(e) for e in site[4::2] ]
    #defect_energies = [ 0.0 for e in site[4::2] ]
    return Site( label, x, [ defect_species[l] for l in defect_labels ], defect_energies, valence=valence )

def format_line( line ):
    """
    Each line in the input file is formatted into separate components to form sites. 
    Args:
        line(str): A line in the input file.
    Returns: 
        line(list): Formatted line from the input file. 
    """
    # x coordinate
    line[1] = float( line[1] )
    line[2] = float( line[2] )
    # defect energies
    for i in range( 4, len(line) ):
        line[i] = float( line[i] )
#        line[i] = 0.0
    return line

def load_site_data( filename, x_min, x_max ):
    """
    Reads in the input data and formats the input data if the x coordinate values are within the calculation region.
    Args:
        filename(string): Filename for importing the input data.
        x_min(float): Minimum x coordinate value defining the calculation region.
        x_max(float): Maximum x coordinate value defining the calculation region. 
   
    Return:
        input_data(list): formatted data for calculation.
    """
    with open( filename, 'r' ) as f:
        input_data = [ line.split() for line in f ]
    input_data = [ format_line( line ) for line in input_data if ( float(line[2]) > x_min and float(line[2]) < x_max ) ]
    return input_data

def mirror_site_data( site_data ):
    """
    Formatted site data is offset so the maximum x coordinate is shifted to an x coordinate of 0.0. The site data with an x coordinate less than 0.0 is mirrored and the shifted and mirrored site data is concatenated to create a system with two grain boundaries. 
    Args:
        site_data(list): Formatted site data.
    Returns:
        site_data_mirrored(list): Formatted site data for a mirrored system. 
    """
    midpoint = max( [ line[1] for line in site_data ] )
    for line in site_data:
        line[1] -= midpoint
    site_data_mirrored = [ copy(l) for l in site_data if l[1] < 0 ]
    for l in site_data_mirrored:
        l[1] = float(l[1]) * -1 
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

def form_continuum_sites( all_sites, x_min, x_max, n_points, b, c, defect_species ):
    """
    NEEDS UPDATING TO WORK WITH NEW CODE FORMAT.
    INCLUDE DOCSTRING
    """
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
