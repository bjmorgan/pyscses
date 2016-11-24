from project.defect_at_site import Defect_at_Site
import numpy as np
from project.constants import fundamental_charge

class Site:
    def __init__( self, x, defect_species, defect_energies, scaling = None ):
        assert( len( defect_species) == len( defect_energies ) )
        self.x = x   
        self.defects = [ Defect_at_Site( d.valence, d.mole_fraction, e, d.fixed ) for d, e in zip( defect_species, defect_energies ) ]
        if scaling:
            self.scaling = scaling
        else:
            self.scaling = np.ones_like( defect_energies )
#       self.defects = [ Defect_Species( valence, mole_fraction ) for  ( valence, mole_fraction ) in defect_data ]
#       self.sites = [ Data(x, energy) for ( x, energy ) in site_data ]

    def sum_of_boltzmann_three( self, phi, temp ):
        return sum( [ d.boltzmann_three( phi, temp ) for d in self.defects ] )

    def probabilities( self, phi, temp ):
        probabilities = []
        for defect in self.defects:
            if defect.fixed:
                probabilities.append( defect.mole_fraction )
            else:  
                probabilities.append( defect.boltzmann_two( phi, temp ) / ( 1.0 + self.sum_of_boltzmann_three( phi, temp ) ) )
        return probabilities
  

    def defect_valences( self ):
        return np.array( [ d.valence for d in self.defects ] )

    def charge( self, phi, temp ):
        '''
        returns the charge in C 
        '''
        #return np.dot( self.probabilities( phi, temp ), self.defect_valences() ) * fundamental_charge
        return np.sum( self.probabilities( phi, temp ) * self.defect_valences() * self.scaling ) * fundamental_charge

    def probabilities_boltz( self, phi, temp ):
        return [ defect.boltzmann_two( phi, temp ) for defect in self.defects ]

    def charge_boltz( self, phi, temp ):
        return np.dot( self.probabilities_boltz( phi, temp ), self.defect_valences() ) * fundamental_charge

    def defect_concentrations( self, phi, temp ):
        return self.density * self.probabilities( phi, temp )
