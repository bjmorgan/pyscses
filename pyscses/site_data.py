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

class SiteData(object):

    def __init__(self,
                 label: str,
                 valence: float,
                 x: float,
                 defect_data: Tuple[DefectData]) -> None:
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

    def as_input_line(self) -> str:
        """Returns a formatted string equivalent to the expected line
        in an input file.

        Args:
            None

        Returns:
            str

        """
        defect_string = " ".join([f"{d.label} {d.energy}" for d in self.defect_data])
        input_line_string = f"{self.label} {self.valence} {self.x} {defect_string}"
        return input_line_string
