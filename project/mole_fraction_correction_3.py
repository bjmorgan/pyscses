from project.defect_species import DefectSpecies
from project.set_of_sites import Set_of_Sites
from project.grid import Grid, delta_x_from_grid
from project.calculation import Calculation
from project.set_up_calculation import calculate_grid_offsets
from bisect import bisect_left
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def mole_fraction_error( input_mole_fractions, target_mole_fractions, temp, mole_fractions, defect_labels, valence, data, bulk_x_min, bulk_x_max, b, c, alpha, conv, boundary_conditions, site_labels ):
    input_mole_fractions = np.array([input_mole_fractions])
    output_mole_fractions = mole_fraction_output( temp, input_mole_fractions, defect_labels, valence, data, bulk_x_min, bulk_x_max, b, c, alpha, conv, boundary_conditions, site_labels )
    o_gd = output_mole_fractions[0,1]
    o_vo = output_mole_fractions[0,0]
    t_gd = target_mole_fractions[0,1]
    t_vo = target_mole_fractions[0,0]
    return ( ( o_gd - t_gd )**2 ) + ( ( o_vo - t_vo )**2 ) 


def mole_fraction_output( temp, mole_fractions, defect_labels, valence, data, bulk_x_min, bulk_x_max, b, c, alpha, conv, boundary_conditions, site_labels ):

    slope_list_vo = []
    intercept_list_vo = []
    slope_list_gd = []
    intercept_list_gd = []

    limits = calculate_grid_offsets( data, bulk_x_min, bulk_x_max )

    for t in temp:
        avg_Vo_molfracs = []
 #       Gd_molfracs = []
        avg_Gd_molfracs = []
   
        for m in mole_fractions:
            Vo_molfracs = []
            Gd_molfracs = []
            print(m, flush=True)
            defect_species = { l : DefectSpecies( l, v, m ) for l, v, m in zip( defect_labels, valence, m) }
            
            all_sites = Set_of_Sites.set_of_sites_from_input_data( data, limits, defect_species ) 
            for site in all_sites.subset( 'site_2' ):
                site.defect_with_label('defect_2').fixed = False
            
            grid = Grid.grid_from_set_of_sites( all_sites, limits[0], limits[1], b, c )
            vo_grid = Grid.grid_from_set_of_sites( all_sites.subset( 'O' ), limits[0], limits[1], b, c )
            gd_grid = grid.grid_from_set_of_sites( all_sites.subset( 'Ce' ), limits[0], limits[1], b, c )
            min_offset = bisect_left(grid.x, bulk_x_min)
            max_offset = bisect_left(grid.x, bulk_x_max) 
            
            calculation_object = Calculation( grid, grid.x[min_offset+2], grid.x[max_offset-2], alpha, conv, t, boundary_conditions )    
            calculation_object.solve()
            calculation_object.form_subgrids(site_labels)
            calculation_object.mole_fractions()
            for mf in calculation_object.mf[site_labels[0]]:
                if mf > 0.0:
                    Vo_molfracs.append(mf)
#            print(Vo_molfracs, calculation_object.mf[site_labels[0]],  len(Vo_molfracs), len(vo_grid.x) )
#            avg_Vo_molfracs.append(np.mean(Vo_molfracs))
            avg_Vo_molfracs.append( (np.sum( Vo_molfracs * delta_x_from_grid( vo_grid.x, limits ) ) / np.sum(delta_x_from_grid( vo_grid.x, limits ) ) ) )
            for mf in calculation_object.mf[site_labels[1]]:
               if mf > 0.0:
                   Gd_molfracs.append(mf)
            avg_Gd_molfracs.append( ( np.sum( Gd_molfracs * delta_x_from_grid( gd_grid.x, limits ) ) / np.sum(delta_x_from_grid( gd_grid.x, limits ) ) ) )
      
    output_mole_fractions = np.array( [ [ a, b ] for a, b in zip( avg_Vo_molfracs, avg_Gd_molfracs ) ] )
    return output_mole_fractions
#    opt_mole_fractions = []

def mole_fraction_correction( temp, mole_fractions, defect_labels, valence, data, bulk_x_min, bulk_x_max, b, c, alpha, conv, boundary_conditions, site_labels ):
    initial_guess = np.array( [ 0.05, 0.2 ] )
    target_mole_fractions = mole_fractions
#    for i, o in zip( mole_fractions, mole_fractions_output ):
    opt_mole_fractions = minimize( mole_fraction_error, initial_guess, args=(  target_mole_fractions, temp, mole_fractions, defect_labels, valence, data, bulk_x_min, bulk_x_max, b, c, alpha, conv, boundary_conditions, site_labels ), bounds = ( (0.0001,1), (0.0001,1) ) )  
    opt_mole_fractions.x = np.array([opt_mole_fractions.x])
    return opt_mole_fractions.x
