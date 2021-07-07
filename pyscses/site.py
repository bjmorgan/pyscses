from __future__ import annotations
import numpy as np
import math
from pyscses.constants import fundamental_charge, boltzmann_eV
from pyscses.grid_point import GridPoint
from pyscses.defect_species import DefectSpecies
from typing import List, Optional, Dict
from pyscses.defect_at_site import DefectAtSite

class LabelError(Exception):
    pass

class Site:
    """The Site class contains all the information about a given site and the defects occupying that site.
    This class contains functions for the calculations which correspond to each individual site, rather than the system as a whole.

    Attributes:
        label (str): Reference label for this site. i.e. 'O' for an oxygen site.
        x (float): x coordinate of the site.
        defect_energies (list): List of segregation energies for all defects present at the site.
        defect_species (list): List of defect species for all defects present at the site.
        defects (list): List of DefectAtSite objects, containing the properties of all individual defects at the site.
        scaling (float): A scaling factor that can be applied in the charge calculation.
        valence (float): The charge of the defect present at the site (in atomic units).
        saturation_parameter (float): Optional saturation parameter as described in
            `Hendricks et al. Sol. Stat. Ionics (2002)`_
            and `Swift et al. Nature Comp. Sci. (2021)`_.
            Setting `saturation_parameter` < `1.0` sets some proportion of excluded sites that are
            unavailable for occupation by any explicit defects.
            Default value is `1.0`, i.e., 100% of sites may be occupied.
        defects (list): List of Defect_Species objects for all defects present at the site.
        sites (list): List containing all x coordinates and corresponding  defect segregation energies.

    .. _Hendricks et al. Sol. Stat. Ionics (2002):
       https://doi.org/10.1016/S0167-2738(02)00484-8

    .. _Swift et al. Nature Comp. Sci. (2021):
       https://doi.org/10.1038/s43588-021-00041-y

    """

    def __init__(self,
                 label: str,
                 x: float,
                 defect_species: List[DefectSpecies],
                 defect_energies: List[float],
                 scaling: Optional[np.ndarray] = None,
                 valence: float = 0.0,
                 saturation_parameter: float = 1.0) -> None:
        """Initialise a Site object.

        Args:
            label (str): Reference label for this site.
            x (float): x coordinate of this site.
            defect_species (list(DefectSpecies)): List of `DefectSpecies` objects (one for each defect species that can occupy this site).
            defect_energies (list(float)): List of defect segregation energies for each defect species at this site.
            scaling (optional, list(float): Optional list of scaling factors for the net charge at this site. Default scaling for each defect species is 1.0.
            valence (optional, float): Optional formal valence for this site in the absence of any defects. Default is 0.0.
            saturation_parameter (optional, float): Optional saturation parameter as described in
                Hendricks et al. Sol. Stat. Ionics (2002) [#HendricksEtAl_SolStatIonics2002]_
                and Swift et al. Nature Comp. Sci. (2021). [#SwiftEtAl_NatureCompSci2021]_.
                A saturation parameter < 1.0 introduces some proportion of excluded sites that are
                unavailable for occupation by any explicit defects. Default is 1.0.

        Raises:
            ValueError if the number of DefectSpecies != the number of defect segregation energies != the number of scaling factors (if passed).

        .. [HendricksEtAl_SolStatIonics2002]:
           https://doi.org/10.1016/S0167-2738(02)00484-8

        .. [SwiftEtAl_NatureCompSci2021]:
           https://doi.org/10.1038/s43588-021-00041-y

        """
        if len(defect_species) != len(defect_energies):
            raise ValueError("len(defect_species) must be equal to len(defect_energies)")
        if scaling:
            if len(defect_species) != len(scaling):
                raise ValueError("len(defect_species) must be equal to len(scaling)")
        self.label = label
        self.x = x
        self.defect_energies = defect_energies
        self.defect_species = defect_species
        self.defects = [DefectAtSite(label=d.label,
                                     valence=d.valence,
                                     mole_fraction=d.mole_fraction,
                                     mobility=d.mobility,
                                     energy=e,
                                     site=self,
                                     fixed=d.fixed)
            for d, e in zip(defect_species, defect_energies)]
        if scaling:
            self.scaling = scaling
        else:
            self.scaling = np.ones_like(defect_energies, dtype=float)
        self.grid_point: Optional[GridPoint] = None
        self.valence = valence
        self.saturation_parameter = saturation_parameter
        self.fixed_defects = tuple(d for d in self.defects if d.fixed)
        self.mobile_defects = tuple(d for d in self.defects if not d.fixed)
        self.alpha = self.saturation_parameter - sum((d.mole_fraction for d in self.fixed_defects))

    def competing_defect_species(self) -> Dict[str, int]:
        """Returns a dictionary reporting the number of fixed and / or mobile defect species that can occupy this site.

        Args:
            None

        Returns
            Dict(str, int): Dictionary {'fixed': n_fixed, 'mobile': n_mobile}

        """
        pass

    def defect_with_label(self,
                          label: str) -> DefectAtSite:
        """Select a defect at this site by the species label.

        Args:
            label (str): Label to identify defect species.

        Returns:
                DefectAtSite: The DefectAtSite that matches the label.

        """
        if not label in (d.label for d in self.defects):
            raise LabelError(f"\"{label}\" does not match any of the defect species labels for this site.")
        else:
            return next(d for d in self.defects if d.label == label)

    def energies(self) -> List[float]:
        """Returns a list of the segregation energies for each defect from self.defects """
        return [d.energy for d in self.defects]

    def average_local_energy(self,
                             method: str = 'mean') -> Optional[np.ndarray]:
        """
        Returns the average local segregation energy for each site based on a specified method.

        Args:
            method (str): The method in which the average segregation energies will be calculated.
                          'mean' - Returns the sum of all values at that site divided by the number of values at that site.
                          'min' - Returns the minimum segregation energy value for that site (appropriate for low temperature calculations).

        Returns:
            numpy.array: Average segregation energies on the site coordinates grid.

        """
        if self.grid_point is not None:
            return self.grid_point.average_site_energy(method)
        else:
            raise ValueError("TODO")

    def probabilities(self,
                      phi: float,
                      temp: float) -> List[float]:
        """Calculates the probabilities of this site being occupied by each defect species.

        Args:
            phi (float): Electrostatic potential at this site in Volts.
            temp (float): Temperature in Kelvin.

        Returns:
            list(float): Probabilities of site occupation on a 1D grid.

        """
        probabilities_dict = {}
        boltzmann_factors = {d.label: d.boltzmann_factor(phi, temp) for d in self.mobile_defects}
        denominator = (self.alpha +
                       sum([d.mole_fraction * (boltzmann_factors[d.label] - 1.0)
                            for d in self.mobile_defects]))
        for defect in self.defects:
            if defect.fixed:
                probabilities_dict[defect.label] = defect.mole_fraction
            else:
                numerator = self.alpha * defect.mole_fraction * boltzmann_factors[defect.label]
                probabilities_dict[defect.label] = numerator / denominator
        return [probabilities_dict[d.label] for d in self.defects]

    def defect_valences(self) -> np.ndarray:
        """Returns an array of valences for each defect from self.defects """
        return np.array([d.valence for d in self.defects])

    def charge(self,
               phi: float,
               temp: float) -> float:
        """
        Charge at this site (in Coulombs).

        Args:
            phi (float):  Electrostatic potential at this site in Volts.
            temp (float): Temperature in Kelvin.

        Returns:
            float: The charge at this site.

        """
        charge =  (self.valence + np.sum(self.probabilities(phi, temp)
                                * self.defect_valences()
                                * self.scaling)) * fundamental_charge
        return charge
