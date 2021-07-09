from copy import copy
from operator import itemgetter
from pyscses.site import Site
import numpy as np
from pyscses.constants import boltzmann_eV
from pyscses.site_data import SiteData, InputFormatError
from bisect import bisect_left, bisect_right
from sklearn.cluster import AgglomerativeClustering # type: ignore
from typing import Tuple, List

def sites_data_from_file(filename: str,
                         x_limits: Tuple[float, float],
                         clustering_threshold: float = 1e-10,
                         site_charge: bool = False) -> List[SiteData]:
    """Reads and pre-processes a set of site data from a file.

    Performs the following operations on the site data:
        1. Excludes any data for sites with x coordinates outside specified limits.
        2. Sorts the site data so that sites are sorted with respect to their x coordinate.
        3. Performs clustering of sites with x coordinates equal within a specified threshold (Default is 0.01 nm).
        4. (optional) Use explicit site charges. If `False` all site charges
            are set to zero. Default is `False`. !!TODO!!

    Args:
        filename (str): The input file.
        x_limits: (tuple(float, float)): x coordinates for the system boundaries.
        clustering_threshold (optional(float)): Distance threshold for clustering sites with similar x coordinated.
            Default is 1e-10.
        site_charge (bool): Set to `True` to use explicit site charges.
            Default is `False.

    Returns:
        list(Site)

    """
    # Read raw data from `filename`:
    with open(filename, 'r') as f:
        input_data = [line.strip() for line in f.readlines() if line]
    # Validate the input data:
    for line_number, line in enumerate(input_data, 1):
        if not SiteData.input_string_is_valid_syntax(line):
            raise InputFormatError(f"Input format error at line number {line_number} in file {filename}")
    # Convert raw data to a list of SiteData objects:
    sites_data = [SiteData.from_input_string(line,
                                             validate_input=False) # We have already validated the input above.
                  for line in input_data]
    # 1. Exclude data for sites with x coordinates outside the specified limits
    sites_data = [sd for sd in sites_data
                  if ((sd.x >= x_limits[0]) and
                      (sd.x < x_limits[1]))]
    # 2. Sort sites data by x coordinate.
    sites_data = sorted(sites_data, key=lambda sd: sd.x)
    # 3. Cluster sites data with similar x coordinates.
    cluster_similar_sites_data(sites_data=sites_data,
                               distance_threshold=1e-10)
    return sites_data

def cluster_similar_sites_data(sites_data: List[SiteData],
                               distance_threshold: float = 1e-10) -> List[SiteData]:
    """Clusters site data with x coordinates equal within a given threshold,
    and adjusts the site data x coordinates to the mean value for each cluster.

    Args:
        sites_data (list(SiteData): List of `SiteData` objects.
        distance_threshold (optional, float): Distance threshold for clustering
            site data with similar x coordinates.
            Default is 1e-10 m.

    Returns:
        (list(SiteData)): Clustered site data.

    """
    if len(sites_data) < 2:
        return sites_data
    all_coordinates = np.array([sd.x for sd in sites_data]).reshape(-1,1)
    cluster = AgglomerativeClustering(n_clusters=None,
                                      affinity='euclidean',
                                      compute_full_tree=True,
                                      linkage='ward',
                                      distance_threshold=distance_threshold)
    cluster.fit(all_coordinates)
    clustered_sites_data = {}
    for label, sd in zip(cluster.labels_, sites_data):
        clustered_sites_data.setdefault(label, []).append(sd)
    for sd_cluster in clustered_sites_data.values():
        mean_x = np.mean([sd.x for sd in sd_cluster])
        for sd in sd_cluster:
            sd.x = mean_x
    return sites_data

def calculate_grid_offsets(filename,
                           x_min,
                           x_max,
                           system):
    """Reads in the input data calculates the distance to the next site outside of the defined calculation region. Allows calculation of the delta_x and volume values for the endmost grid points.

    Args:
        filename(string): Filename for importing the input data.
        x_min(float): Minimum x coordinate value defining the calculation region.
        x_max(float): Maximum x coordinate value defining the calculation region.
        system(str): 'single' for single grain boundary systems, 'double' for systems where the input data has been mirrored to give two grain boundaries.

    Returns:
        list, list: distance between the midpoint of the endmost sites and the midpoint of the next site outside of the calculation region for the first and last sites respectively. Distance between the endmost sites and the next site outside of the calculation region for the first and last sites respectively.

    """
    with open( filename, 'r' ) as f:
        input_data = [ line.split() for line in f ]
        input_data = [ format_line( line, 0.0, 0.0 ) for line in input_data ]
        input_data = sorted( input_data, key=lambda x: x[2] )
        input_data = cluster_similar_sites( input_data )
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

