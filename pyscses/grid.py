from __future__ import annotations
import numpy as np
import math
from pyscses.constants import boltzmann_eV
from bisect import bisect_left
from scipy.interpolate import griddata # type: ignore
from typing import Union, Tuple, List
from pyscses.grid_point import GridPoint
from pyscses.defect_species import DefectSpecies
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pyscses.set_of_sites import SetOfSites

def phi_at_x(phi: np.ndarray,
             coordinates: np.ndarray,
             x: float) -> float:
    """
    Assigns each site x coordinate a grid point and returns the electrostatic potential at the grid point closest to the x coordinate.

    Args:
        phi (np.array): electrostatic potential on 1D grid.
        coordinates (np.array): 1D grid of ordered numbers over a region.
        x (float): Site x coordinate.

    Returns:
        float: The electrostatic potential at the x coordinate with position [index].
    """
    index = index_of_grid_at_x(coordinates, x)
    return phi[index]

def energy_at_x(energy: np.ndarray,
                coordinates: np.ndarray,
                x: float) -> float:
    """
    Assigns each site x coordinate a grid point and returns the segregation energy at the grid point closest to the x coordinate.

    Args:
        energy (np.array): Segregation energies on 1D grid.
        coordinates (np.array): 1D grid of ordered numbers over a region.
        x (float): Site x coordinate.

    Returns:
        energy[index] (float): The segregation energy at the x coordinate with position [index].
    """

    index = index_of_grid_at_x(coordinates, x)
    return energy[index]

def index_of_grid_at_x(coordinates: np.ndarray,
                       x: float) -> int:
    """
    Assigns each site x coordinate to a position on a regularly or irregularly spaced grid.
    Returns the index of the grid point clostest to the value x

    Args:
        coordinates (np.array): 1D grid of ordered numbers over a region.
        x (float): Site x coordinate

    Returns:
        int: Index of grid position closest to the site x coordinate.
    """
    return closest_index(coordinates, x)

def closest_index(myList: Union[list[float], np.ndarray],
                  myNumber: float) -> int:
    """
    Assumes myList is sorted. Returns index of closest value to myNumber.
    If two numbers are equally close, return the index of the smallest number.

    Args:
        myList  (list): List of numbers to compare against. This should be sorted.
        myNumber (float): The number to compare against myList.

    Returns:
        pos (int): Index of position of number in myList which is closest to myNumber.
    """
    myList = list(myList)
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return (len(myList) - 1)
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos - 1

def delta_x_from_grid(coordinates: np.ndarray,
                      limits: Tuple[float, float]) -> np.ndarray:
    """
    Calculates the distance between the midpoints of each consecutive site.
    Inserts the calculated distance to the next grid point outside of the calculation region to
    the first and last position as the delta_x value for the endmost sites.

    Args:
        coordinates (np.array): 1D grid of ordered numbers over a region.
        limits (list): Distance between the midpoint of the endmost sites and the midpoint of the respective site outside of the calculation region.

    Returns:
        delta_x (np.array): Distance between the midpoints of each consecutive site.
    """
    delta_x = np.zeros_like(coordinates)
    delta_x = (coordinates[2:] - coordinates[:-2]) / 2.0
    delta_x = np.insert(delta_x, 0, limits[0])
    delta_x = np.insert(delta_x, len(delta_x), limits[1])
    return delta_x

class Grid:

    def __init__(self,
                 x_coordinates: np.ndarray,
                 b: float,
                 c: float,
                 limits: Tuple[float, float],
                 limits_for_laplacian: Tuple[float, float],
                 set_of_sites: SetOfSites) -> None:
       # x_coordinates need to be sorted for the delta_x calculation in volumes_from_grid
        """
        Grid objects contain information and methods for the calculation grid.

        Args:
            delta_x (np.array): Distance between the midpoint of each consecutive grid point.
            volumes (np.array): Volume of each consecutive grid point.
            points (list): GridPoint object at each grid point in Grid.
            x (np.array): The x-coordinates for each grid point.
            limits(list): distance between the midpoint of the endmost sites and the midpoint of the next site outside of the calculation region for the first and last sites respectively.
            limits_for_laplacian(list): distance between the endmost sites and the next site outside of the calculation region for the first and last sites respectively.
            b (float):                b dimension for every grid-point.
            c (float):                c dimension for every grid-point.
            limits ([float,float])    x-coordinates for the minimum and maximum grid edges.
            set_of_sites (SetOfSites): Set of Site objects that populate the Grid.
            defect_Species (list): Defect species that populate the Grid.
        Returns:
            None
        """
        self.delta_x = delta_x_from_grid(x_coordinates, limits)
        self.volumes = self.delta_x * b * c
        self.x = x_coordinates
        self.points = [GridPoint(x=x, volume=v)
                       for x, v in zip(self.x, self.volumes)]
        self.limits = limits
        self.limits_for_laplacian = limits_for_laplacian
        self.b = b
        self.c = c
        self.set_of_sites = set_of_sites
        self.defect_species: List[DefectSpecies] = []
        for site in set_of_sites:
            i = index_of_grid_at_x(self.x, site.x)
            self.points[i].sites.append(site)
            site.grid_point = self.points[i]
            for defect_species in site.defect_species:
                if defect_species not in self.defect_species:
                    self.defect_species.append(defect_species)

    def __getitem__(self,
                    key: int) -> GridPoint:
        return self.points[key]

    def charge(self,
               phi: np.ndarray,
               temp: float) -> np.ndarray:
        """
        Calculates the overall charge at each point on a grid.

        Args:
            phi (np.array): Electrostatic potential on a 1D grid.
            temp (float): Temperature in Kelvin.

        Returns:
            np.array: Overall charge at each point on a 1D grid.
        """
        charge = np.zeros_like(self.x)
        for i, point in enumerate(self.points):
            charge[i] = sum([site.charge(
                                phi=phi_at_x(
                                    phi=phi,
                                    coordinates=self.x,
                                    x=site.x),
                                temp=temp)
                             for site in point.sites])
        return np.array(charge)

    def rho(self,
            phi: np.ndarray,
            temp: float) -> np.ndarray:
        """
        Calculates charge density at each point on a grid, by dividing the overall charge at that grid point by the grid point volume.

        Args:
            phi (np.array): Electrostatic potential on a 1D grid.
            temp (float): Temperature in K.

        Returns:
            np.array: The charge density at each point on a grid.
        """
        rho = np.zeros_like(self.x)
        rho = self.charge(phi, temp) / self.volumes
        return rho

# BEN Is this ever called?
#     def defect_valences(self) -> np.ndarray:
#         """
#             (np.array): Valences for each defect from self.defects
#         Returns:
#         for site in self.points.sites:
#         """
#             return np.array([d.valence for d in site.defects])


#     def average_site_energies(self,
#                               method: str = 'mean') -> np.ndarray:
#         """
#         Returns the average segregation energy for all grid points based on a specified method
#
#         Args:
#             method (str): The method in which the average segregation energies will be calculated.
#                           'mean' - Returns the sum of all values at that site divided by the number of values at that site.
#                           'min' - Returns the minimum segregation energy value for that site (appropriate for low temperature calculations).
#
#         Returns:
#             np.array: Average segregation energies on a 1D grid.
#         """
#         return np.array([p.average_site_energy(method) for p in self.points]).T

#     def interpolated_energies(self) -> List[float]:
#         """
#         Returns:
#            list: The average site energies linearly interpolated onto a regularly spaced grid
#         """
#
#         energies = []
#         for row in self.average_site_energies():
#             x_new, e_new = np.array([(x, e) for x, e in zip(self.x, row) if e]).T
#             interpolated_energies = griddata(x_new, e_new, self.x, fill_value = 0.0)
#             energies.append(interpolated_energies)
#         return energies

    def subgrid(self,
                subset_species: str) -> Grid:
        """
        Creates a sub-grid for a single defect species.

        Args:
            subset_species (str): Site species to separate into sub-grid.

        Returns:
            :obj:`Grid`: Grid object for subset of data.
        """
        return Grid.from_set_of_sites(set_of_sites=self.set_of_sites.subset(subset_species),
                                      limits=self.limits,
                                      limits_for_laplacian=self.limits_for_laplacian,
                                      b=self.b,
                                      c=self.c)

    @classmethod
    def from_set_of_sites(cls: object,
                          set_of_sites: SetOfSites,
                          limits: Tuple[float, float],
                          limits_for_laplacian: Tuple[float, float],
                          b: float,
                          c: float) -> Grid:
        """
        Creates a grid from a given Set_of_Sites object.

        Args:
            set_of_sites (SetOfSites): Set_of_Sites object containing a set of all Site objects.
            limits (list): distance between the midpoint of the endmost sites and the midpoint of the next site outside of the calculation region for the first and last sites respectively.
            limits_for_laplacian (list): distance between the endmost sites and the next site outside of the calculation region for the first and last sites respectively.
            b (float):                b dimension for every grid-point.
            c (float):                c dimension for every grid-point.

        Returns:
            :obj:`Grid`: Grid object for the given set of sites.

        """
        return Grid(x_coordinates=np.unique([site.x for site in set_of_sites]),
                    b=b,
                    c=c,
                    limits=limits,
                    limits_for_laplacian=limits_for_laplacian,
                    set_of_sites=set_of_sites)
