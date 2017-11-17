from project.defect_species import DefectSpecies
from project.set_of_sites import Set_of_Sites
from project.grid import Grid, delta_x_from_grid
from project.calculation import Calculation
from project.set_up_calculation import calculate_grid_offsets
from bisect import bisect_left
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def mole_fraction_correction( temp, mole_fractions, defect_labels, valence, data, bulk_x_min, bulk_x_max, b, c, alpha, conv, boundary_conditions, site_labels ):

    slope_list_vo = []
    intercept_list_vo = []
    slope_list_gd = []
    intercept_list_gd = []

    limits = calculate_grid_offsets( data, bulk_x_min, bulk_x_max )

    for t in temp:
        avg_Vo_molfracs = []
        avg_Gd_molfracs = []
        for m in mole_fractions:
            Vo_molfracs = []

            defect_species = { l : DefectSpecies( l, v, m ) for l, v, m in zip( defect_labels, valence, m) }
            
            all_sites = Set_of_Sites.set_of_sites_from_input_data( data, limits, defect_species ) 
            for site in all_sites.subset( 'site_2' ):
                site.defect_with_label('defect_2').fixed = False
            
            grid = Grid.grid_from_set_of_sites( all_sites, limits[0], limits[1], b, c )
            vo_grid = Grid.grid_from_set_of_sites( all_sites.subset( 'O' ), limits[0], limits[1], b, c )
     
            min_offset = bisect_left(grid.x, bulk_x_min)
            max_offset = bisect_left(grid.x, bulk_x_max) 
            
            calculation_object = Calculation( grid, grid.x[min_offset+2], grid.x[max_offset-2], alpha, conv, t, boundary_conditions )    
            calculation_object.solve()
            calculation_object.form_subgrids(site_labels)
            calculation_object.mole_fractions()
            for mf in calculation_object.mf[site_labels[0]]:
                if mf > 0.0:
                    Vo_molfracs.append(mf)
            avg_Vo_molfracs.append( ( np.sum( Vo_molfracs * delta_x_from_grid( vo_grid.x, limits ) ) / np.sum(delta_x_from_grid( vo_grid.x, limits ) ) ) )
        slope_vo, intercept_vo, rvalue, pvalue, stderr = stats.linregress( mole_fractions[:,0], avg_Vo_molfracs )
        input_vo_mf = []
        for m in mole_fractions[:,0]:
            input_vo_mf.append( ( m - intercept_vo ) / slope_vo )

        mole_fractions = np.array([[a, m] for a, m in zip(input_vo_mf, mole_fractions[:,1])])
        for m in mole_fractions:
            Gd_molfracs = []

            defect_species = { l : DefectSpecies( l, v, m ) for l, v, m in zip( defect_labels, valence, m) }
            
            all_sites = Set_of_Sites.set_of_sites_from_input_data( data, limits, defect_species ) 
            for site in all_sites.subset( 'site_2' ):
                site.defect_with_label('defect_2').fixed = False
            
            grid = Grid.grid_from_set_of_sites( all_sites, limits[0], limits[1], b, c )
            gd_grid = Grid.grid_from_set_of_sites( all_sites.subset( 'Ce' ), limits[0], limits[1], b, c )
     
            min_offset = bisect_left(grid.x, bulk_x_min)
            max_offset = bisect_left(grid.x, bulk_x_max) 
            
            calculation_object = Calculation( grid, grid.x[min_offset+2], grid.x[max_offset-2], alpha, conv, t, boundary_conditions )    
            calculation_object.solve()
            calculation_object.form_subgrids(site_labels)
            calculation_object.mole_fractions()
            for mf in calculation_object.mf[site_labels[1]]:
                if mf > 0.0:
                     Gd_molfracs.append(mf)
            avg_Gd_molfracs.append( ( np.sum( Gd_molfracs * delta_x_from_grid( gd_grid.x, limits ) ) / np.sum(delta_x_from_grid( gd_grid.x, limits ) ) ) )
        slope_gd, intercept_gd, rvalue, pvalue, stderr = stats.linregress( mole_fractions[:,1], avg_Gd_molfracs )
        input_gd_mf = []
        for m in mole_fractions[:,1]:
            input_gd_mf.append( ( m - intercept_gd ) / slope_gd )
        mole_fractions = np.array( [ [m, a] for m, a in zip(mole_fractions[:,0], input_gd_mf) ] )
        return mole_fractions





 
