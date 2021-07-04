from __future__ import annotations
import math
from pyscses.constants import boltzmann_eV
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pyscses.site import Site

class DefectAtSite:

    """
    The DefectAtSite class contains the information about each defect at each site, its valence, mole fraction and the segregation energy for that defect at that site.

    The methods in this class combine to give the correct statistics for site occupation in solid electrolytes, derived from the electrochemical potentials of a Fermi-Dirac like distribution. This term has been split into three functions for simplicity. The resulting equations take into account that the probability that a site is occupied by a defect can not exceed 1, and also accounts for competition of like charged defects.
    A Mott-Schottky approximation can be enforced. The defects can be fixed to their bulk mole fractions throughout the system, equivalent to assuming that the defects are immobile and not allowed to redistribute through the system.

    Attributes:

        label (string): string describing the defect species i.e. 'Vo' for an oxygen vacancy.
        valence (float): The charge of the defect, in atomic units.
        mole_fraction (float): The bulk mole fraction of the defect present in the system.
        mobility (float): The bulk mobility of the defect species. Default = 0.0. TODO Where is this default mobility set?
        energy (float): The segregation energy for the defect when occupying the respective site.
        site (:obj:`Site`): The pyscses.Site object that corresponding to each defect at each x coordinate.
        fixed (bool): set whether this defect species can redistribute to an equilibrium distribution. Default=False.

    """
    def __init__(self,
                 label: str,
                 valence: float,
                 mole_fraction: float,
                 mobility: float,
                 energy: float,
                 site: Site,
                 fixed: bool = False) -> None:
        self.label = label
        self.valence = valence
        self.mole_fraction = mole_fraction
        self.mobility = mobility
        self.energy = energy
        self.site = site
        self.fixed = fixed

    def potential_energy(self,
                         phi: float) -> float:
        """
        Potential energy for the defect at this site.

        Args:
            phi (float): electrostatic potential at this site.

        Returns:
            float: The electrochemical potential.

        """
        return (phi * self.valence) + self.energy

    def boltzmann_one(self,
                      phi: float,
                      temp: float) -> float:
        r"""Boltzmann statistics calculation - part one

        .. math:: \exp\left(\frac{-\Phi z+\Delta E}{k_BT}\right)

        Args:
            phi (float): Electrostatic potential.
            temp (float): Temperature in Kelvin.

        Returns:
            float: Boltzmann factor.

        """
        return math.exp(-self.potential_energy(phi) / (boltzmann_eV * temp))

    def boltzmann_two(self,
                      phi: float,
                      temp: float) -> float:
        r"""Boltzmann statistics calculation - part two

        .. math:: x\exp\left(\frac{-\Phi z+\Delta E}{K_BT}\right)

        Args:
            phi (float): Electrostatic potential.
            temp (float): Temperature of calculation.

        Returns:
            float: Boltzmann statistics * mole fraction

        """
        return self.mole_fraction * self.boltzmann_one(phi, temp)

    def boltzmann_three(self,
                        phi: float,
                        temp: float) -> float:
        r"""Boltzmann statistics calculation - part three

        .. math:: x\left(\exp\left(-\frac{\Phi z+\Delta E}{K_BT}\right)-1\right)

        Args:
            phi (float): Electrostatic potential.
            temp (float): Temperature of calculation.

        Returns:
            float: ( Boltzmann statistics - 1 ) * mole fraction.

        """
        return self.mole_fraction * (self.boltzmann_one(phi, temp) - 1.0)
