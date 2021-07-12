from __future__ import annotations

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
