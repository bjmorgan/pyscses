from project.defect_at_site import Defect_at_Site
import numpy as np
import math
from project.constants import fundamental_charge, boltzmann_eV

class Site:
    """ The site class contains all the information about a given site and the defect occupying that site.
        This class contains functions for the calculations which correspond to each individual site, rather than the whole system. """
    def __init__( self, label, x, defect_species, defect_energies, scaling = None, valence = 0 ):
        assert( len( defect_species) == len( defect_energies ) )
        self.label = label
        self.x = x   
        self.defect_species = defect_species
        self.defects = [ Defect_at_Site( d.label, d.valence, d.mole_fraction, e, self, d.fixed ) for d, e in zip( defect_species, defect_energies ) ]
        if scaling:
            self.scaling = scaling
        else:
            self.scaling = np.ones_like( defect_energies )
        self.grid_point = None
        self.valence = valence
#       self.defects = [ Defect_Species( valence, mole_fraction ) for  ( valence, mole_fraction ) in defect_data ]
#       self.sites = [ Data(x, energy) for ( x, energy ) in site_data ]

    def defect_with_label( self, label ):
        """ Returns a list of defects which correspond to the given label """
        return [ d for d in self.defects if d.label == label ][0]

    def energies( self ):
        """ Returns a list of the segregation energies for each defect from self.defects """
        return [ d.energy for d in self.defects ]

    def average_local_energy( self, method = 'mean' ):
        """ 
        Returns the average local segregation energy for each site based on a specified method 

        Args: 
            method (str): The method in which the average segregation energies will be calculated.
                          'mean' - Returns the sum of all values at that site divided by the number of values at that site.
                          'min' - Returns the minimum segregation energy value for that site (appropriate for low temperature calculations).

        Returns:
            average site energies (np.array): Average segregation energies on the site coordinates grid.
  
        """
        return self.grid_point.average_site_energy( method )

    def sum_of_boltzmann_three( self, phi, temp ):
        """
        Calculates the sum of the calculated boltzmann_three values for each defect at each site.  

        Args: 
            phi (float): Electrostatic potential at the site.
            temp (float): Temperature of calculation in Kelvin.

        Returns:
            The sum of ( mole fraction * exp ^ -( ( phi * z )/kBT ) - 1 ) for each defect at each site. 
        
        """
        return sum( [ d.boltzmann_three( phi, temp ) for d in self.defects ] )

    def probabilities( self, phi, temp ):
        """
        Calculates the probability of each site being occupied by a given defect using Fermi-Dirac like statistics.

        Args:
            phi (floaty):   Electrostatic potential at this site.
            temp (float):   Temperature of calculation.

        Returns: 
            probabilities (list): Probabilities of site occupation on a 1D grid. 
        """
        probabilities = []
        for defect in self.defects:
            if defect.fixed:
                probabilities.append( defect.mole_fraction )
            else:  
                probabilities.append( defect.boltzmann_two( phi, temp ) / ( 1.0 + self.sum_of_boltzmann_three( phi, temp ) ) )
        return probabilities 

    def defect_valences( self ):
        """ Returns an array of valences for each defect from self.defects """
        return np.array( [ d.valence for d in self.defects ] )


    def charge( self, phi, temp ):
        """
        Calculates the charge in Coulombs at each site under Fermi-Dirac like statistics.

        Args:
            phi (float):  Electrostatic potential at this site.
            temp (float): Temperature of calculation.

        Returns:
            charge (np.array): The charge on a 1D grid.
        """
        charge =  ( self.valence +  np.sum( self.probabilities( phi, temp ) * self.defect_valences() * self.scaling ) ) * fundamental_charge
        return charge

    def probabilities_boltz( self, phi, temp ):
        """
    
        Calculates the probability of each site being occupied by a given defect using Boltzmann statistics.

        Args:
            phi (float): Electrostatic potential at this site.
            temp (float): Temperature of calculation.

        Returns: 
            boltzmann_probabilities (list): Probabilities of site occupation on a 1D grid. 

        """
        boltzmann_probabilities = [ defect.boltzmann_two( phi, temp ) for defect in self.defects ]
        return boltzmann_probabilities

    def charge_boltz( self, phi, temp ):
        """
        Calculates the charge in Coulombs at each site under Boltzmann statistics.

        Args:
            phi (float): Electrostatic potential at this site
            temp (float): Temperature of calculation.

        Returns:
            charge (np.array): The charge on a 1D grid.
        """
        charge =  np.sum( self.probabilities_boltz( phi, temp ) * self.defect_valences() * self.scaling ) * fundamental_charge
        return charge


