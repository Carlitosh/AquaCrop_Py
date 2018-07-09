#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from hydrology_funs import *

import logging
logger = logging.getLogger(__name__)

class Drainage(object):
    """Class to infiltrate incoming water"""
    
    def __init__(self, Drainage_variable):
        self.var = Drainage_variable

    def initial(self):
        self.var.FluxOut = np.zeros((self.var.nComp, self.var.nRotation, self.var.nLat, self.var.nLon))
        self.var.DeepPerc = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        
    def dynamic(self):
        """Function to redistribute stored soil water"""
        th, DeepPerc, FluxOut = drainage(self.var.th,
                                         self.var.th_s_comp,
                                         self.var.th_fc_comp,
                                         self.var.th_fc_adj,
                                         self.var.ksat_comp,
                                         self.var.tau_comp,
                                         self.var.dz,
                                         self.var.dzsum)
        self.var.th  = th
        self.var.DeepPerc = DeepPerc
        self.var.FluxOut = FluxOut
