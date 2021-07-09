from copy import copy
from operator import itemgetter
from pyscses.site import Site
import numpy as np
from pyscses.constants import boltzmann_eV
from pyscses.site_data import SiteData, InputFormatError
from bisect import bisect_left, bisect_right
from sklearn.cluster import AgglomerativeClustering # type: ignore
from typing import Tuple, List, Dict

def sites_data_from_file(filename: str,
                         x_limits: Tuple[float, float],
                         clustering_threshold: float = 1e-10,
                         site_charge: bool = False) -> Tuple[List[SiteData], Tuple[SiteData, SiteData]]:
    """Reads and pre-processes a set of site data from a file.

    Performs the following operations on the site data:
        1. Sorts the site data so that sites are sorted with respect to their x coordinate.
        2. Performs clustering of sites with x coordinates equal within a specified threshold (Default is 0.01 nm).
        3. (optional) Use explicit site charges. If `False` all site charges
            are set to zero. Default is `False`. !!TODO!!

    Args:
        filename (str): The input file.
        x_limits: (tuple(float, float)): x coordinates for the system boundaries.
        clustering_threshold (optional(float)): Distance threshold for clustering sites with similar x coordinated.
            Default is 1e-10.
        site_charge (bool): Set to `True` to use explicit site charges.
            Default is `False.

    Returns:
        list(SiteData), tuple(SiteData, SiteData): TODO Explain what is being returned.

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
    # 1. Sort sites data by x coordinate.
    sites_data = sorted(sites_data, key=lambda sd: sd.x)
    # 2. Cluster sites data with similar x coordinates.
    cluster_similar_sites_data(sites_data=sites_data,
                               distance_threshold=1e-10)
    # 3. (Optional). Set site charges to zero if site_charge == True.
    # TODO
    # 4. Select sites with x coordinates within the specified x limits,
    #    plus the two adjacent sites outside the upper and lower bounds.
    x_coords = [sd.x for sd in sites_data]
    index_lower = np.searchsorted(x_coords, x_limits[0])
    index_upper = np.searchsorted(x_coords, x_limits[1])
    bounding_sites_data = sites_data[index_lower-1], sites_data[index_upper]
    sites_data = sites_data[index_lower:index_upper]
    return sites_data, bounding_sites_data

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
    clustered_sites_data: Dict[int, List[SiteData]] = {}
    for label, sd in zip(cluster.labels_, sites_data):
        clustered_sites_data.setdefault(label, []).append(sd)
    for sd_cluster in clustered_sites_data.values():
        mean_x = np.mean([sd.x for sd in sd_cluster])
        for sd in sd_cluster:
            sd.x = mean_x
    return sites_data
