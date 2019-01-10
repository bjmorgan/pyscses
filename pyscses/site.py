from pyscses.defect_at_site import Defect_at_Site
import numpy as np
import math
from pyscses.constants import fundamental_charge, boltzmann_eV

class Site:
    """ The site class contains all the information about a given site and the defect occupying that site.
        This class contains functions for the calculations which correspond to each individual site, rather than the system as a whole.

    Args:
        label (str): refers to what the defect is called i.e. 'Vo' for an oxygen vacancy. 
        x (float): x coordinate of the site. 
	defect_energies (list): List of segregation energies for all defects present at the site.
        defect_species (list): List of defect species for all defects present at the site.
	defects (list): List of Defect_at_Site objects, containing the properties of all individual defects at the site.
	scaling (float): A scaling factor that can be applied in the charge calculation.
	valence (float): The charge of the defect present at the site (in atomic units).
	defects (list): List of Defect_Species objects for all defects present at the site.
	sites (list): List containing all x coordinates and corresponding  defect segregation energies.
    """

    def __init__( self, label, x, defect_species, defect_energies, scaling = None, valence = 0 ):
        assert( len( defect_species) == len( defect_energies ) )
        self.label = label
        self.x = x  
        self.defect_energies = defect_energies 
        self.defect_species = defect_species
        self.defects = [ Defect_at_Site( d.label, d.valence, d.mole_fraction, d.mobility, e, self, d.fixed ) for d, e in zip( defect_species, defect_energies ) ]
        if scaling:
            self.scaling = scaling
        else:
            self.scaling = np.ones_like( defect_energies )
        self.grid_point = None
        self.valence = valence
#       self.defects = [ Defect_Species( valence, mole_fraction ) for  ( valence, mole_fraction ) in defect_data ]
#       self.sites = [ Data(x, energy) for ( x, energy ) in site_data ]

    def defect_with_label( self, label ):
        """
	Returns a list of defects which correspond to the given label
	
	Args:
	    label (str): Label to identify defect species.
    
	Returns:
	    (list): List of Defect_at_Site objects for a specific defect species. 
	"""
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
        Calculates the sum of the calculated boltzmann_three values for each defect at each site.  i

	.. math:: \sum(x_i\exp(-\Phi_xz/kT)-1)

        Args: 
            phi (float): Electrostatic potential at the site.
            temp (float): Temperature of calculation in Kelvin.

        Returns:
            (float): The sum of Boltzmann terms.
        
        """
        return sum( [ d.boltzmann_three( phi, temp ) for d in self.defects ] )

    def probabilities( self, phi, temp ):
        """
        Calculates the probability of each site being occupied. Derived from the chemical potential term for a Fermi-Dirac like distribution.

        Args:
            phi (float):   Electrostatic potential at this site.
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
        Calculates the overall charge in Coulombs at each site.

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
    
        Calculates the probability of each site being occupied by a given defect. Derived from the chemical potential including a Boltzmann distribution.

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
        Calculates the charge in Coulombs at each site using Boltzmann probabilities.

        Args:
            phi (float): Electrostatic potential at this site
            temp (float): Temperature of calculation.

        Returns:
            charge (np.array): The charge on a 1D grid.
        """
        charge =  np.sum( self.probabilities_boltz( phi, temp ) * self.defect_valences() * self.scaling ) * fundamental_charge
        return charge


