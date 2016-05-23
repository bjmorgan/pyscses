import math
from project.constants import boltzmann_eV

class Defect_at_Site:
    def __init__( self, valence, mole_fraction, energy, fixed = False ):
        self.valence = valence
        self.mole_fraction = mole_fraction
        self.energy = energy
        self.fixed = fixed 

    def boltzmann_one( self, phi, temp ):
        return math.exp( -( ( phi * self.valence ) + self.energy ) / ( boltzmann_eV * temp ) )

    def boltzmann_two( self, phi, temp ):
        return self.mole_fraction * self.boltzmann_one( phi, temp )

    def boltzmann_three( self, phi, temp ):
        return self.mole_fraction * ( self.boltzmann_one( phi, temp ) - 1.0 )
