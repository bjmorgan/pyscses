from project.defect_species import Defect_Species
from project.set_of_sites import Set_of_Sites
from project.grid import Grid
from project.calculation import Calculation

import numpy as np
from scipy import stats

def mole_fraction_correction( temp, mole_fractions, defect_labels, valence, data, limits, b, c, alpha, conv, boundary_conditions, site_labels ):

    slope_list = []
    intercept_list = []

    for t in temp:
        Vo_molfracs = []
        avg_Vo_molfracs = []
    
        for m in mole_fractions:
        
            defect_species = { l : Defect_Species( l, v, m ) for l, v, m in zip( defect_labels, valence, m) }
    
            all_sites = Set_of_Sites.set_of_sites_from_input_data( data, limits, defect_species )
            for site in all_sites.subset( 'site_2' ):
                site.defect_with_label('defect_2').fixed = True
        
            grid = Grid.grid_from_set_of_sites( all_sites, limits[0], limits[1], b, c )
       
    
            calculation_object = Calculation( grid, alpha, conv, t, boundary_conditions )
            calculation_object.solve()
            calculation_object.form_subgrids(site_labels)
            calculation_object.mole_fractions
            for mf in calculation_object.mf[site_labels[0]]:
                if mf > 0.0:
                    Vo_molfracs.append(mf)
            avg_Vo_molfracs.append(np.mean(Vo_molfracs))
        slope, intercept, rvalue, pvalue, stderr = stats.linregress( mole_fractions[:,0], avg_Vo_molfracs )
        slope_list.append( slope )
        intercept_list.append( intercept )
    
    percentage_Gd = 20
    desired_mobile_defect_MF = ( percentage_Gd / 100 ) / 4
    mole_fractions = np.array( [ MF( desired_mobile_defect_MF, s, i ) for s, i in zip( slope_list, intercept_list ) ] )

    return mole_fractions
