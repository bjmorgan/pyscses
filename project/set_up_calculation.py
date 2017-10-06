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
    x = float(site[2])
    defect_labels = site[3::2]
    defect_energies = [ float(e) for e in site[4::2] ]
#    defect_energies = [ 0.0 for e in site[4::2] ]
    return Site( label, x, [ defect_species[l] for l in defect_labels ], defect_energies, valence=valence )

def format_line( line ):
    """
    another docstring
    """
    # x coordinate
    line[1] = float( line[1] )
    line[2] = float( line[2] )
    # defect energies
    for i in range( 4, len(line) ):
        line[i] = float( line[i] )
    return line

def load_site_data( filename, x_min, x_max ):
    """
    docstring goes here.
    """
    with open( filename, 'r' ) as f:
        input_data = [ line.split() for line in f ]
    input_data = [ format_line( line ) for line in input_data if ( float(line[2]) > x_min and float(line[2]) < x_max ) ]
    return input_data

def mirror_site_data( site_data ):
    """
    docstring
    """
    midpoint = max( [ line[1] for line in site_data ] )
    for line in site_data:
        line[1] -= midpoint
    site_data_mirrored = [ copy(l) for l in site_data if l[1] < 0 ]
    for l in site_data_mirrored:
        l[1] = float(l[1]) * -1 
    return site_data + site_data_mirrored

#def new_load_site_data( filename, x_min, x_max ):
#    """
#    """
#    with open( filename, 'r' ) as f:
#        input_data = [ line.split() for line in f ]
#        x_coords = np.array( [ l[2] for l in input_data ] )
#        min_index = bisect_right( x_coords, x_min )
#        max_index = bisect_right( x_coords, x_max )
#        return [ format_line( line ) for line in input_data[ min_index: max_index ] ]

def calculate_grid_offsets( filename, x_min, x_max ):
    with open( filename, 'r' ) as f:
        input_data = [ line.split() for line in f ]
        x_coords = np.unique( np.array( [ float(l[2]) for l in input_data ] ) )
        min_index = bisect_left( x_coords, x_min )
#        min_offset = x_coords[min_index] - x_coords[min_index-1]
        max_index = bisect_right( x_coords, x_max )
#        max_offset = x_coords[max_index] - x_coords[max_index-1]
        limits = [ x_coords[min_index-1], x_coords[max_index] ]
#    return min_offset, max_offset
    return limits
