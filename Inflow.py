#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

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
        self.GwIn = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        
    def dynamic(self):
        """Function to calculate capillary rise in the presence of a 
        shallow groundwater table
        """

        # Initialize groudwater inflow array
        GwIn = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

        # Water table in soil profile: calculate horizontal inflow; get
        # groundwater table elevation on current day
        zBot = np.cumsum(self.var.dz_xy, axis=0)
        zTop = zBot - self.var.dz_xy
        zMid = (zTop + zBot) / 2

        # For compartments below water table, set to saturation
        dth = np.zeros((self.var.nCrop, self.var.nComp, self.var.nLat, self.var.nLon))
        cond1 = (np.broadcast_to(self.var.WTinSoil[:,None,:,:], self.var.th.shape)
                 & (zMid >= np.broadcast_to(self.var.zGW[:,None,:,:], self.var.th.shape)))
        cond11 = (cond1 & (self.var.th < self.var.th_s_comp))

        # Update water content
        dth[cond11] = (self.var.th_s_comp - self.var.th)[cond11]
        self.var.th[cond11] = self.var.th_s_comp[cond11]

        # Update groundwater inflow
        GwIn_comp = dth * 1000 * self.var.dz_xy
        self.var.GwIn = np.sum(GwIn_comp, axis=0)
