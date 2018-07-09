#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from hydrology_funs import *

import logging
logger = logging.getLogger(__name__)

class Inflow(object):
    """Class to check for presence of a groundwater table and, if 
    present, to adjust compartment water contents and field 
    capacities where necessary
    """
    def __init__(self, GroundwaterInflow_variable):
        self.var = GroundwaterInflow_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.GwIn = np.copy(arr_zeros)
        
    def dynamic(self):
        """Function to calculate capillary rise in the presence of a 
        shallow groundwater table
        """
        self.var.th, self.var.GwIn = groundwater_inflow(
            self.var.th, self.var.th_s_comp,
            self.var.WTinSoil,
            self.var.zGW,
            self.var.dz,
            self.var.dzsum)
            
        # print self.var.th[0,0,0,0]
