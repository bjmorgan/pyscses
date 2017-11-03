from project.defect_species import DefectSpecies
from project.set_of_sites import Set_of_Sites
from project.grid import Grid
from project.calculation import Calculation
from project.set_up_calculation import calculate_grid_offsets

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def MF( desired_mobile_defect_mf, slope, intercept ):
    """
    Calculates the mole fraction that is required in the input to get the desired output mole fraction.
    Due to some numerical noise, when the simulation is run the defect mole fractions vary from what they should be and this code uses a solpe and intercept from linear regression of the input and output mole fractions to find what the input mole fraction should be to achieve the desired output.

    Args:
        desired_mobile_defect_MF (float): desired output mole fraction for the mobile defect.
        slope (float): slope calculated from linear regression.
        intercept (float): intercept calcuated from linear regression.

    Returns:
        MF (float): input mole fraction to achieve desired output mole fraction.
    """
    mobile_defect_mf = ( desired_mobile_defect_mf - intercept ) / slope
    MF = [ (mobile_defect_mf), ( mobile_defect_mf * 4 ) ]
    return MF

def mole_fraction_correction( temp, mole_fractions, defect_labels, valence, data, bulk_x_min, bulk_x_max, b, c, alpha, conv, boundary_conditions, site_labels ):

    slope_list = []
    intercept_list = []

    limits = calculate_grid_offsets( data, bulk_x_min, bulk_x_max )
    limits_2 = calculate_grid_offsets( data, limits[0], limits[1] )

    for t in temp:
        Vo_molfracs = []
        avg_Vo_molfracs = []
    
        for m in mole_fractions:
            print(m, flush=True)
            defect_species = { l : DefectSpecies( l, v, m ) for l, v, m in zip( defect_labels, valence, m) }
    
            all_sites = Set_of_Sites.set_of_sites_from_input_data( data, limits, defect_species )
            for site in all_sites.subset( 'site_2' ):
                site.defect_with_label('defect_2').fixed = False
        
            grid = Grid.grid_from_set_of_sites( all_sites, limits_2[0], limits_2[1], b, c )
       
    
            calculation_object = Calculation( grid, bulk_x_min, bulk_x_max, alpha, conv, t, boundary_conditions )
            calculation_object.solve()
            calculation_object.form_subgrids(site_labels)
            calculation_object.mole_fractions()
            for mf in calculation_object.mf[site_labels[0]]:
                if mf > 0.0:
                    Vo_molfracs.append(mf)
            avg_Vo_molfracs.append(np.mean(Vo_molfracs))
        slope, intercept, rvalue, pvalue, stderr = stats.linregress( mole_fractions[:,0], avg_Vo_molfracs )
        slope_list.append( slope )
        intercept_list.append( intercept )
        plt.plot(mole_fractions[:,0], avg_Vo_molfracs)    
    percentage_Gd = 20
    desired_mobile_defect_MF = ( percentage_Gd / 100 ) / 4
    mole_fractions = np.array( [ MF( desired_mobile_defect_MF, s, i ) for s, i in zip( slope_list, intercept_list ) ] )

    return mole_fractions

#def phi_correction( t, mole_fractions, defect_labels, valence, data, limits, b, c, alpha, conv, boundary_conditions, site_labels ):
#    defect_species = { l : DefectSpecies( l, v, m ) for l, v, m in zip( defect_labels, valence, mole_fractions) }

#    all_sites = Set_of_Sites.set_of_sites_from_input_data( data, limits, defect_species )
#    for site in all_sites.subset( 'site_2' ):
#        site.defect_with_label('defect_2').fixed = False

#    grid = Grid.grid_from_set_of_sites( all_sites, limits[0], limits[1], b, c )


#    calculation_object = Calculation( grid, alpha, conv, t, boundary_conditions, 0.0 )
#    calculation_object.solve()
#    return ( np.sum( calculation_object.phi * delta_x_from_grid( grid.x, grid.limits ) ) / np.sum(delta_x_from_grid( grid.x, grid.limits ) ) )
