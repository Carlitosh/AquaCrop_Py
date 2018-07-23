#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class Germination(object):

    def __init__(self, Germination_variable):
        self.var = Germination_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.DelayedGDDs = np.copy(arr_zeros)
        self.var.DelayedCDs  = np.copy(arr_zeros)
        self.var.Germination = np.copy(arr_zeros.astype(bool))
        self.var.AgeDays     = np.copy(arr_zeros)
        self.var.AgeDays_NS  = np.copy(arr_zeros)

    def reset_initial_conditions(self):
        self.var.DelayedGDDs[self.var.GrowingSeasonDayOne] = 0
        self.var.DelayedCDs[self.var.GrowingSeasonDayOne] = 0
        self.var.Germination[self.var.GrowingSeasonDayOne] = False
        self.var.AgeDays[self.var.GrowingSeasonDayOne] = 0
        self.var.AgeDays_NS[self.var.GrowingSeasonDayOne] = 0
        
    def dynamic(self):
        """Function to check if crop has germinated"""
        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()
            
        # Here we force zGerm to have a maximum value equal to the depth of the
        # deepest soil compartment
        zgerm = np.copy(self.var.zGerm)
        zgerm[zgerm > np.sum(self.var.dz, axis=0)] = np.sum(self.var.dz, axis=0)

        # Find compartments covered by top soil layer affecting germination
        comp_sto = (np.round(self.var.dzsum_xy * 1000) <= np.round(np.broadcast_to(zgerm[:,None,:,:], self.var.th.shape) * 1000))  # round to nearest mm

        # Calculate water content in top soil layer
        arr_zeros = np.zeros((self.var.nCrop, self.var.nComp, self.var.nLat, self.var.nLon))
        Wr_comp   = np.copy(arr_zeros)
        WrFC_comp = np.copy(arr_zeros)
        WrWP_comp = np.copy(arr_zeros)

        # Determine fraction of compartment covered by top soil layer
        factor = 1. - np.round(((self.var.dzsum_xy - zgerm) / self.var.dz_xy), 3)
        factor = np.clip(factor, 0, 1) * comp_sto

        # Increment water storages (mm)
        Wr_comp = np.round((factor * 1000 * self.var.th * self.var.dz_xy))
        Wr_comp = np.clip(Wr_comp, 0, None)
        Wr = np.sum(Wr_comp, axis=1)

        WrFC_comp = np.round((factor * 1000 * self.var.th_fc_comp * self.var.dz_xy))
        WrFC = np.sum(WrFC_comp, axis=1)

        WrWP_comp = np.round((factor * 1000 * self.var.th_wp_comp * self.var.dz_xy))
        WrWP = np.sum(WrWP_comp, axis=1)

        # Calculate proportional water content
        WrTAW = WrFC - WrWP
        WcProp = 1 - np.divide((WrFC - Wr), WrTAW, out=np.zeros_like(WrTAW), where=WrTAW!=0)
        
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
        cond7 = (self.var.DAP > self.var.MaxCanopyCD) & self.var.GrowingSeasonIndex  # NB not originally in this function
        self.var.AgeDays_NS[cond7] = (self.var.DAP - self.var.MaxCanopyCD)[cond7]
        
        self.var.Germination[np.logical_not(self.var.GrowingSeasonIndex)] = False
        self.var.DelayedCDs[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.DelayedGDDs[np.logical_not(self.var.GrowingSeasonIndex)] = 0
