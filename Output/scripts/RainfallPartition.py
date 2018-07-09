#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from hydrology_funs import *

import logging
logger = logging.getLogger(__name__)

class RainfallPartition(object):
    """Class to infiltrate incoming water"""
    
    def __init__(self, RainfallPartition_variable):
        self.var = RainfallPartition_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLon, self.var.nLat))        
        self.Runoff = np.copy(arr_zeros)
        self.Infl = np.copy(arr_zeros)
        
    def dynamic(self):
        """Function to partition rainfall into surface runoff and 
        infiltration using the curve number approach.
        """
        Runoff, Infl = rainfall_partition(self.var.precipitation,
                                          self.var.th,
                                          self.var.th_fc_comp,
                                          self.var.th_wp_comp,
                                          self.var.zCN,
                                          self.var.AdjCN,
                                          self.var.CN,
                                          self.var.CNbot,
                                          self.var.CNtop,
                                          self.var.dz,
                                          self.var.dzsum,
                                          self.var.Bunds,
                                          self.var.zBund)
        self.var.Runoff = Runoff
        self.var.Infl = Infl

