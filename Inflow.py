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
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.GwIn = np.copy(arr_zeros)
        
    def dynamic(self):
        """Function to calculate capillary rise in the presence of a 
        shallow groundwater table
        """
        # self.var.th, self.var.GwIn = groundwater_inflow(
        #     self.var.th, self.var.th_s_comp,
        #     self.var.WTinSoil,
        #     self.var.zGW,
        #     self.var.dz,
        #     self.var.dzsum)
        dz = self.var.dz[:,None,None,None] * np.ones((self.var.nRotation, self.var.nLat, self.var.nLon))[None,:,:,:]
        dzsum = self.var.dzsum[:,None,None,None] * np.ones((self.var.nRotation, self.var.nLat, self.var.nLon))[None,:,:,:]

        # Initialize groudwater inflow array
        GwIn = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))

        # Water table in soil profile: calculate horizontal inflow; get
        # groundwater table elevation on current day
        zBot = np.cumsum(dz, axis=0)
        zTop = zBot - dz
        zMid = (zTop + zBot) / 2

        # For compartments below water table, set to saturation
        dth = np.zeros((self.var.nComp, self.var.nRotation, self.var.nLat, self.var.nLon))
        cond1 = (self.var.WTinSoil & (zMid >= self.var.zGW))
        cond11 = (cond1 & (self.var.th < self.var.th_s_comp))

        # Update water content
        dth[cond11] = (self.var.th_s_comp - self.var.th)[cond11]
        self.var.th[cond11] = self.var.th_s_comp[cond11]

        # Update groundwater inflow
        GwIn_comp = dth * 1000 * dz
        self.var.GwIn = np.sum(GwIn_comp, axis=0)
