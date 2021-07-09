import numpy as np
from pyscses.set_up_calculation import sites_data_from_file
from typing import Tuple

class InputData(object):
	
    def __init__(self,
                 filename: str,
                 x_limits: Tuple[float, float],
                 system: str):
        self.sites_data, self.adjacent_sites_data = sites_data_from_file(filename=filename,
                                                                         x_limits=x_limits)
        self.x_limits = x_limits
        self.site_x_coords = np.unique([sd.x for sd in self.sites_data])
        self.system = system

    def limits(self) -> Tuple[float, float]:
        min_offset = (self.site_x_coords[1] - self.adjacent_sites_data[0].x)/2.0
        max_offset = (self.adjacent_sites_data[1].x - self.site_x_coords[-2])/2.0
        if self.system == 'single':
            return (min_offset, max_offset)
        elif self.system == 'double':
            return (min_offset, min_offset)
        else:
            raise Exception()

    def limits_for_laplacian(self) -> Tuple[float, float]:
        min_offset = (self.site_x_coords[1] - self.adjacent_sites_data[0].x)/2.0
        max_offset = (self.adjacent_sites_data[1].x - self.site_x_coords[-2])/2.0
        if self.system == 'single':
            return (self.site_x_coords[0] - self.adjacent_sites_data[0].x,
                    self.adjacent_sites_data[1].x - self.site_x_coords[-1])
        elif self.system == 'double':
                return (self.site_x_coords[0] - self.adjacent_sites_data[0].x,
                        self.site_x_coords[0] - self.adjacent_sites_data[0].x)
        else:
            raise Exception()
