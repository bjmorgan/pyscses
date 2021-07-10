from __future__ import annotations
from collections import namedtuple
from typing import Tuple
import re
from pyscses.defect_data import DefectData

class InputFormatError(Exception):
    pass

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
        self.defect_labels = set(dd.label for dd in self.defect_data)

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
    def from_input_string(cls,
                          input_string: str,
                          validate_input: bool = True) -> SiteData:
        """Parse a formatted string in the input file format and return
        a corresponding `SiteData` instance.

        Args:
            input_string (str): String describing the input data for this site.
            validate_input (optional, bool): Flag for whether to validate the input string.
                Default is `True`.

        Return:
            SiteData

        Raises:
            ValueError: if the input is found to be invalid syntax.

        """
        if validate_input:
            if not cls.input_string_is_valid_syntax(input_string):
                raise InputFormatError(f"Invalid syntax for site data: {input_string}")
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

    @staticmethod
    def input_string_is_valid_syntax(string: str) -> bool:
        """Test whether a given input string is a valid format to construct a `SiteData` object.

        Args:
            string (str): The string to be tested.

        Returns
            bool

        """
        input_re = re.compile(("(\w+)"                     # label (string)
                               "\s+"                       # whitespace
                               "([+-\.\d]+)"               # valence (float)
                               "\s+"                       # whitespace
                               "([+-\.e\d]+)"              # x-coordinate (float)
                               "(\s+(\w+)\s+([+-\.e\d]+))+" # 1 or more [defect_label defect_energy] pairs
                                                           # (string, float)
                               "\Z"))                      # end of string
        return bool(input_re.match(string))
