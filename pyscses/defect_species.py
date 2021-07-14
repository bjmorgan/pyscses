class DefectSpecies(object):
    """
    The DefectSpecies class describes the properties for a single defect species present in the system.

    Attributes:
        label (string): Label for this defect species e.g. "Vo" for an oxygen vacancy.
        valence (float): The formal charge for this defect species, in atomic units.
        mole_fraction (float): The bulk mole fraction of this defect species.
        mobility (float): Mobility of this defect species. Default is `0.0`.
        can_equilibrate (bool): Specifies whether this defect species is allowed to redistribute to achieve an equilibrium distribution,
        or is kept at its input distribution. Default is `True`.

    """

    def __init__(self,
                 label: str,
                 valence: float,
                 mole_fraction: float,
                 mobility: float = 0.0,
                 can_equilibrate: bool = True) -> None:
        if not isinstance(label, str):
            raise TypeError("When initialising a DefectSpecies object, the label argument must be a string.")
        if not isinstance(valence, float):
            raise TypeError("When initialising a DefectSpecies object, the valence argument must be a float.")
        if not isinstance(mole_fraction, float):
            raise TypeError("When initialising a DefectSpecies object, the mole_fraction argument must be a float.")
        self.label = label
        self.valence = valence
        self.mole_fraction = mole_fraction
        self.mobility = mobility
        self.can_equilibrate = can_equilibrate
