class DefectSpecies:
    """ 
    The DefectSpecies class includes the information about each defect species present in the system. 

    Attributes:

        label (string): refers to what the defect is called i.e. 'Vo' for an oxygen vacancy. 
        valence (float): The charge of the defect, in atomic units.
        mole_fraction (float): The mole fraction of the defect present in the system. 
        fixed (bool): set whether this defect species can redistribute to an equilibrium distriution. Default=False.
        mobility (float): mobility of this defect species. Default = 0.0.
    """

    def __init__( self, label, valence, mole_fraction, fixed=False, mobility=0.0 ):
        assert( isinstance( label, str ) )
        assert( isinstance( valence, float ) )
        assert( isinstance( mole_fraction, float ) )
        self.label = label
        self.valence = valence
        self.mole_fraction = mole_fraction
        self.fixed = fixed
        self.mobility = mobility
