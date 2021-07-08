from __future__ import annotations
from collections import namedtuple
from typing import Tuple

class DefectData(object):

    def __init__(self,
                 label: str,
                 energy: float) -> None:
        """Instantiate a DefectData object

        Args:
            label (str): Label for this defect species.
            energy (float): Segregation energy for this defect.

        Returns:
            None

        """
        self.label = label
        self.energy = energy

    def __eq__(self,
               other: object) -> bool:
        if not isinstance(other, DefectData):
            return NotImplemented
        return (self.label == other.label) and (self.energy == other.energy)

class SiteData(object):

    def __init__(self,
                 label: str,
                 valence: float,
                 x: float,
                 defect_data: Tuple[DefectData, ...]) -> None:
        """Instantiate a SiteData object.

        Args:
            label (str): Site label used to reference this site.
            valence (float): Site formal charge.
            x (float): Site x coordinate.
            defect_data (tuple(DefectData)): Data describing each defect species that can occupy this site.

        Returns:
            None

        """
        self.label = label
        self.valence = valence
        self.x = x
        self.defect_data = defect_data

    def as_input_string(self) -> str:
        """Returns a formatted string equivalent to the expected line
        in an input file.

        Args:
            None

        Returns:
            str

        """
        defect_string = " ".join([f"{d.label} {d.energy}" for d in self.defect_data])
        input_string = f"{self.label} {self.valence} {self.x} {defect_string}"
        return input_string

    @classmethod
    def from_input_string(self,
                          input_string: str) -> SiteData:
        """Parse a formatted string in the input file format and return
        a corresponding `SiteData` instance.

        Args:
            input_string (str): String describing the input data for this site.

        Return:
            SiteData

        """
        input = input_string.split()
        label = input[0]
        valence = float(input[1])
        x = float(input[2])
        defect_labels = input[3::2]
        defect_energies = [float(s) for s in input[4::2]]
        defect_data = tuple(DefectData(label=l, energy=e)
                            for l, e in zip(defect_labels, defect_energies))
        site_data = SiteData(label=label,
                             valence=valence,
                             x=x,
                             defect_data=defect_data)
        return site_data
