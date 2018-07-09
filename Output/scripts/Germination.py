#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from crop_growth_funs import *

import logging
logger = logging.getLogger(__name__)

class Germination(object):

    def __init__(self, Germination_variable):
        self.var = Germination_variable

    def initial(self):
        pass

    def water_content_affecting_germination(self):

        dims = self.var.th_fc_comp.shape
        nc, nr, nlat, nlon = dims[0], dims[1], dims[2], dims[3]

        # Here we force zGerm to have a maximum value equal to the depth of the
        # deepest soil compartment
        zgerm = np.copy(self.var.zGerm)
        zgerm[zgerm > np.sum(self.var.dz, axis=0)] = np.sum(dz, axis=0)

        # Add rotation, lat, lon dimensions to dz and dzsum
        dz = self.var.dz[:,None,None,None] * np.ones((nr, nlat, nlon))
        dzsum = self.var.dzsum[:,None,None,None] * np.ones((nr, nlat, nlon))

        # Find compartments covered by top soil layer affecting germination
        comp_sto = (np.round(dzsum * 1000) <= np.round(zgerm * 1000))  # round to nearest mm

        # Calculate water content in top soil layer
        arr_zeros = np.zeros((nc, nr, nlat, nlon))
        Wr_comp   = np.copy(arr_zeros)
        WrFC_comp = np.copy(arr_zeros)
        WrWP_comp = np.copy(arr_zeros)

        # Determine fraction of compartment covered by top soil layer
        factor = 1. - np.round(((dzsum - zgerm) / dz), 3)
        factor = np.clip(factor, 0, 1) * comp_sto

        # Increment water storages (mm)
        Wr_comp = np.round((factor * 1000 * self.var.th * dz))
        Wr_comp = np.clip(Wr_comp, 0, None)
        Wr = np.sum(Wr_comp, axis=0)

        WrFC_comp = np.round((factor * 1000 * self.var.th_fc_comp * dz))
        WrFC = np.sum(WrFC_comp, axis=0)

        WrWP_comp = np.round((factor * 1000 * self.var.th_wp_comp * dz))
        WrWP = np.sum(WrWP_comp, axis=0)

        # Calculate proportional water content
        WrTAW = WrFC - WrWP
        WcProp = 1 - np.divide((WrFC - Wr), WrTAW, out=np.zeros_like(WrTAW), where=WrTAW!=0)
        return WcProp
    
    def dynamic(self):
        """Function to check if crop has germinated"""
        WcProp = water_content_affecting_germination(
            self.var.th,
            self.var.th_fc_comp,
            self.var.th_wp_comp,
            self.var.dz,
            self.var.dzsum,
            self.var.zGerm)
        
        # Check if water content is above germination threshold
        cond4 = (self.var.GrowingSeasonIndex & (WcProp >= self.var.GermThr) & (np.logical_not(self.var.Germination)))
        self.var.Germination[cond4] = True

        # Increment delayed growth time counters if germination is yet to occur
        cond5 = (self.var.GrowingSeasonIndex & (np.logical_not(self.var.Germination)))
        self.var.DelayedCDs[cond5] += 1
        self.var.DelayedGDDs[cond5] += self.var.GDD[cond5]

        # Update ageing days counter
        DAPadj = (self.var.DAP - self.var.DelayedCDs)
        cond6 = (DAPadj > self.var.MaxCanopyCD) & self.var.GrowingSeasonIndex
        self.var.AgeDays[cond6] = (DAPadj - self.var.MaxCanopyCD)[cond6]
        
        self.var.Germination[np.logical_not(self.var.GrowingSeasonIndex)] = False
        self.var.DelayedCDs[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.DelayedGDDs[np.logical_not(self.var.GrowingSeasonIndex)] = 0
