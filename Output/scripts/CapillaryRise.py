#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from hydrology_funs import *

import logging
logger = logging.getLogger(__name__)

class CapillaryRise(object):
    def __init__(self, CapillaryRise_variable):
        self.var = CapillaryRise_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLon, self.var.nLat))
        self.var.CrTot  = np.copy(arr_zeros)
        
    def dynamic(self):
        """Function to calculate capillary rise from a shallow 
        groundwater table
        """
        if self.var.WaterTable:

            # Get maximum capillary rise for bottom compartment
            zBot = np.sum(self.var.dz)              
            zBotMid = zBot - (self.var.dz[-1] / 2)  # depth to midpoint of bottom layer
            MaxCR = maximum_capillary_rise(
                self.var.ksat[-1,:], self.var.aCR[-1,:],
                self.var.bCR[-1,:], self.var.zGW, zBotMid)

            # Check for restrictions on upward flow caused by properties of
            # compartments that are not modelled in the soil water balance

            # Find top of next soil layer that is not within modelled soil profile
            zTopLayer = np.cumsum(self.var.zLayer[np.arange(0, (self.var.layerIndex[-1] + 1))])
            layeri = self.var.layerIndex[-1]  # layer number of bottom compartment
            LimCR = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))

            # TODO: should this be nlayer, rather than nlayer-1
            while np.any(zTopLayer < self.var.zGW) & (layeri < (self.var.nLayer - 1)):
                layeri += 1
                LimCR = maximum_capillary_rise(
                    self.var.ksat[layeri,:], self.var.aCR[layeri,:],
                    self.var.bCR[layeri,:], self.var.zGW, zTopLayer)
                MaxCR = np.clip(MaxCR, None, LimCR)
                zTopLayer += self.var.zLayer[layeri]

            thnew, CrTot = capillary_rise_fun(
                self.var.th, self.var.th_fc_comp, self.var.th_fc_adj,
                self.var.th_wp_comp, self.var.fshape_cr_comp,
                self.var.aCR_comp, self.var.bCR_comp, MaxCR,
                self.var.FluxOut, self.var.zGW, self.var.dz)
            
            self.var.th = np.copy(thnew)
            self.var.CrTot = CrTot

        # print self.var.th[0,0,0,0]
