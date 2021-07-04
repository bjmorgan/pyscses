from __future__ import annotations
from scipy.interpolate import griddata # type: ignore
import numpy as np
import math
from pyscses.set_up_calculation import site_from_input_file, load_site_data
from pyscses.grid import index_of_grid_at_x, phi_at_x, energy_at_x
from pyscses.constants import boltzmann_eV
from pyscses.defect_species import DefectSpecies
from bisect import bisect_left
from typing import List, Iterator, Tuple, Optional
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pyscses.site import Site
    from pyscses.grid import Grid

class SetOfSites:
    """The SetOfSites object groups together all of the Site objects into one object and
    contains functions for the calculations that provide properties of all of the sites
    together rather than individually.

    Attributes:
        sites (list(Site)): List of Site objects.

    """
    def __init__(self,
                 sites: List[Site]) -> None:
        """Initialise a SetOfSites object."""
        self.sites = sites

    def __add__(self,
                other: SetOfSites) -> SetOfSites:
        """Allows the concatenation of multiple SetOfSites objects"""
        if type(other) is not SetOfSites:
            raise TypeError
        return SetOfSites(self.sites + other.sites)

    def __getitem__(self,
                    index: int) -> Site:
        """Returns the site corresponding to a given index """
        return self.sites[index]

    def __iter__(self) -> Iterator[Site]:
        """Iterator over self.sites"""
        return iter(self.sites)

    def subset(self,
               label: str) -> SetOfSites:
        """Returns a subset of all the sites which contain a particular defect """
        return SetOfSites(sites=[s for s in self.sites if s.label == label])

    def get_coords(self,
                   label: str) -> List[float]:
        """
        Returns a list of the x coordinates for all the sites wich contain a particular defect

    	Args:
    	    label (str): Label identifying the required defect species.

    	Returns:
    	    list(float): List of site coordinates for a specific defect species.

    	"""
        return [s.x for s in self.sites if s.label == label]

# BEN: Is this used?
#     def calculate_energies_on_grid(self,
#                                    grid: Grid,
#                                    phi: np.ndarray) -> np.ndarray:
#         """
#         Returns an array of energies at their points on a one-dimensional grid.
#
#         Args:
#             grid (Grid): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
#             phi (array): electrostatic potential on a one-dimensional grid.
#
#         Returns:
#             array: energies at their grid points
#
#         """
#         energies_on_grid = np.zeros_like(grid.x)
#         for site in self.sites:
#             energies_on_grid[index_of_grid_at_x(grid.x, site.x)] =+ energy_at_x(site.defect_energies, grid.x, site.x)
#         return energies_on_grid

    def calculate_probabilities(self,
                                grid: Grid,
                                phi: np.ndarray,
                                temp: float) -> np.ndarray:
        """
        Calculates the probability of a site being occupied by its corresponding defect.

        Args:
            grid (Grid): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            phi (array): electrostatic potential on a one-dimensional grid.
            temp (float): Absolute temperature.

        Returns:
            array: probabilities of defects occupying each site using their grid points

        """
        probability = np.zeros_like(grid.x)
        for i,j in enumerate(grid.x):
            prob = []
            for site in self.sites:
                if j == site.x:
                    prob.append(site.probabilities(phi_at_x( phi, grid.x, site.x ), temp))
            if len(prob) == 0:
                probability[i] = 0
            else:
                probability[i] = np.mean(prob)
        return probability

    def calculate_defect_density(self,
                                 grid: Grid,
                                 phi: np.ndarray,
                                 temp: float) -> np.ndarray:
        """
        Calculates the defect density at each site.

        Args:
            grid (Grid): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates.
            phi (array): electrostatic potential on a one-dimensional grid.
            temp (float): Absolute temperature.

        Returns:
            array: defect density for each site using their grid points

        """
        defect_density = np.zeros_like( grid.x )
        for site in self.sites:
            i = index_of_grid_at_x( grid.x, site.x )
            defect_density[ i ] += np.asarray( site.probabilities( phi_at_x( phi, grid.x, site.x ), temp ) ) / grid.volumes[ i ]
        return defect_density

    def subgrid_calculate_defect_density(self,
                                         sub_grid: Grid,
                                         full_grid: Grid,
                                         phi: np.ndarray,
                                         temp: float) -> np.ndarray:
        """
        Calculates the defect density at each site for a given subset of sites.

        Args:
            subgrid (Grid): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates. For a given subset of sites.
            full_grid (Grid): Grid object - contains properties of the grid including the x coordinates and the volumes. Used to access the x coordinates. For all sites.
            phi (array): electrostatic potential on a one-dimensional grid.
            temp (float): Absolute temperature.

        Returns:
            array: defect density for each site using their grid points

        """
        defect_density = np.zeros_like( sub_grid.x )
        for site in self.sites:
            i = index_of_grid_at_x( sub_grid.x, site.x )
            defect_density[ i ] += np.asarray( site.probabilities( phi_at_x( phi, full_grid.x, site.x ), temp ) ) / sub_grid.volumes[ i ]
        return defect_density

# BEN: Does not appear to be called?
#     def form_continuum_sites(all_sites,
#                              x_min,
#                              x_max,
#                              n_points,
#                              b,
#                              c,
#                              defect_species,
#                              limits_for_laplacian,
#                              site_labels,
#                              defect_labels):
#         """
#         Creates a SetOfSites object for sites interpolated onto a regular grid, this is equivalent to assuming a continuum approximation.
#
#         Args:
#             all_sites (SetOfSites): Orginal SetOfSites object from full data.
#             x_min (float): Minimum x coordinate value defining the calculation region.
#             x_max (float): Maximum x coordinate value defining the calculation region.
#             n_points (int): Number of points that the data should be interpolated on to.
#             b (float): b dimension for every grid point.
#             c (float): c dimension for every grid point.
#             defect_species (object): Class object containing information about the defect species present in the system.
#             limits for laplacian (list): distance between the endmost sites and the midpoint of the next site outside of the calculation region for the first and last sites respectively.
#             site_labels( list ): List of strings for the different site species.
#             defect_labels (list): List of strings for the different defect species.
#
#         Returns:
#             :obj:`SetOfSites`: Sites interpolated onto a regular grid.
#
# 	"""
#
#         grid = np.linspace( x_min, x_max, n_points )
#         limits = [grid[1] - grid[0], grid[1] - grid[0]]
#         sites = []
#         for label, d_label in zip(site_labels, defect_labels):
#             scaling = len( all_sites.subset( label ) ) / len( grid )
#             continuum_grid = Grid(grid, b, c, limits, limits_for_laplacian, all_sites.subset(label))
#             average_energies = np.array( [ site.average_local_energy( method = 'mean' )[0] for site in all_sites.subset( label ) ] )
#             new_energies = griddata( ( [ site.x for site in all_sites.subset( label ) ] ), average_energies, grid, method = 'nearest' )
#             for x, e in zip( grid, new_energies):
#                 sites.append( Site( label, x, [ defect_species[ d_label ] ], [e], scaling = np.array( scaling ) ) )
#         return SetOfSites( sites ), limits

    @ classmethod
    def set_of_sites_from_input_data(cls: object,
                                     filename: str,
                                     limits: Tuple[float, float],
                                     defect_species: List[DefectSpecies],
                                     site_charge: bool,
                                     core: str,
                                     temperature: float,
                                     offset: float = 0.0) -> SetOfSites:
        """
        Takes the data from the input file and creates a SetOfSites object for those sites.
	    The input data file is a .txt file where each line in the file corresponds to a site. The values in each line are formatted and separated into the corresponding properties before creating a Site object for each site.

        Args:
    	    filename (str): Name of the input file to be parsed.
    	    limits (list): Minimum and maximum x coordinated defining the calculation region.
                defect_species (object): Class object containing information about the defect species present in the system.
                site_charge (bool): The site charge refers to the contribution to the overall charge of a site given by the original, non-defective species present at that site. True if the site charge contribution is to be included in the calculation, False if it is not to be included.
    	    core (str): Core definition. 'single' = Single segregation energy used to define the core. 'multi-site' = Layered segregation energies used to define the core while the energies fall in the region of positive and negative kT. 'all' = All sites between a minimum and maximum x coordinate used in calculation.
    	    temperature (float): Temperature that the calculation is being run at.

    	Returns:
    	    :obj:`SetOfSites`: `SetOfSites` object for the input data.

     	"""
        site_data = load_site_data(filename, limits[0], limits[1], site_charge, offset)
        energies = [line[4] for line in site_data]
        min_energy = min(energies)
        if core == 'single':
            for line in site_data:
                if line[4] > min_energy:
                    line[4] = 0.0
        if core == 'multi_site':
            for line in site_data:
                if ( -boltzmann_eV * temperature) <= line[4] <= ( boltzmann_eV * temperature ):
                    line[4] = 0.0
        return SetOfSites([site_from_input_file(line, defect_species, site_charge, core, temperature) for line in site_data])

    @ classmethod
    def core_width_analysis( cls, input_data, limits, defect_species, site_charge, core, temperature ):
        """
	Calculated the width of the 'core' region. This is given as the region where the segregation energies in the system are within a region of positive to negative kT.

        Args:
	    input_data (file): A .txt file where each line includes information about a site.
	    limits (list): Minimum and maximum x coordinates defining the calculation region.
	    defect_species (object): Class object containing information about the defect species present in the system.
            site_charge (bool): The site charge refers to the contribution to the overall charge of a site given by the original, non-defective species present at that site. True if the site charge contribution is to be included in the calculation, False if it is not to be included.
	    core (str): Core definition. Allowed keywords: 'single' = Single segregation energy used to define the core. 'multi-site' = Layered segregation energies used to define the core while the energies fall in the region of positive and negative kT. 'all' = All sites between a minimum and maximum x coordinate used in calculation.
	    temperature (float): Temperature that the calculation is being run at.

	Returns:
	    float: Distance between the minimum and maximum x coordinates where the segregation energy is in the range of positive to negative kT.

	"""
        site_data = load_site_data( input_data, limits[0], limits[1], site_charge )
        energies = [ line[4] for line in site_data ]
        min_energy = min(energies)
        if core == 'single':
            for line in site_data:
                if line[4] > min_energy:
                    line[4] = 0.0
        #print(boltzmann_eV * temperature, flush=True)
        if core == 'multi_site':
            for line in site_data:
                if ( -boltzmann_eV * temperature) <= line[4] <= ( boltzmann_eV * temperature ):
                    line[4] = 0.0
        energies = [line[4] for line in site_data ]
        x = [line[2] for line in site_data ]
        x_seg = np.column_stack(( x, energies ))
        minval = np.min(x_seg[:,0][np.nonzero(x_seg[:,1])])
        maxval = np.max(x_seg[:,0][np.nonzero(x_seg[:,1])])
        core_width = maxval-minval
        return core_width
