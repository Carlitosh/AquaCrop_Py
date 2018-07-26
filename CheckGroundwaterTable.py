#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class CheckGroundwaterTable(object):
    """Class to check for presence of a groundwater table and, if 
    present, to adjust compartment water contents and field 
    capacities where necessary
    """
    def __init__(self, CheckGroundwaterTable_variable):
        self.var = CheckGroundwaterTable_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.th_fc_adj = np.copy(self.var.th_fc_comp)
        self.var.WTinSoil = np.copy(arr_zeros.astype(bool))

    def reset_initial_conditions(self):
        self.var.WTinSoil[self.var.GrowingSeasonDayOne] = False
        
    def dynamic(self):

        # reset initial conditions
        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()
        
        if self.var.WaterTable:
            # Copy depth to groundwater, and add crop dimension for convenience
            self.var.zGW = self.var.zGW[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
            zGW_comp = self.var.zGW[:,None,:,:] * np.ones((self.var.nComp))[None,:,None,None]

            # get the mid point of each compartment
            zBot = np.cumsum(self.var.dz)
            zTop = zBot - self.var.dz
            zMid = (zTop + zBot) / 2
            zMid = zMid[None,:,None,None] * np.ones((self.var.nCrop,self.var.nLat,self.var.nLon))[:,None,:,:]

            # Check if water table is within modelled soil profile
            WTinSoilComp = (zMid >= zGW_comp)
            # WTinSoilComp = (zMid >= self.var.zGW)
            self.var.th[WTinSoilComp] = self.var.th_s_comp[WTinSoilComp]

            # Flatten WTinSoilComp to provide an array with dimensions
            # (ncrop, nLat, nLon), indicating crops where the water
            # table is in the soil profile
            self.var.WTinSoil = np.sum(WTinSoilComp, axis=1).astype(bool)

            # get Xmax
            Xmax = np.zeros((self.var.nCrop,self.var.nComp,self.var.nLat,self.var.nLon))
            cond1 = self.var.th_fc_comp <= 0.1
            cond2 = self.var.th_fc_comp >= 0.3
            cond3 = np.logical_not(cond1 | cond2) # i.e. 0.1 < fc < 0.3
            Xmax[cond1] = 1
            Xmax[cond2] = 2
            pF = 2 + 0.3 * (self.var.th_fc_comp - 0.1) / 0.2
            Xmax_cond3 = np.exp(pF * np.log(10)) / 100
            Xmax[cond3] = Xmax_cond3[cond3]
            cond4 = (zGW_comp < 0) | ((zGW_comp - zMid) >= Xmax)
            # cond4 = (self.var.zGW < 0) | ((self.var.zGW - zMid) >= Xmax)

            # Index of the compartment to which each element belongs (shallow ->
            # deep, i.e. 1 is the shallowest)
            compartment = (np.arange(1, self.var.nComp + 1)[None,:,None,None] * np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))[:,None,:,:])

            # Index of the lowest compartment (i.e. the maximum value) for which
            # cond4 is met, cast to all compartments (achieved by multiplying
            # compartments by cond4 to set elements that do not equal the
            # condition to zero, but retain the compartment number of elements
            # that do meet the condition
            cond4_max_compartment = (np.amax(compartment * cond4, axis=1)[:,None,:,:] * np.ones((self.var.nComp))[None,:,None,None])

            # Now, identify compartments that are shallower than the deepest
            # compartment for which cond4 is met
            cond4 = (compartment <= cond4_max_compartment)

            # 'cond4' is a special case because if ANY compartment meets the
            # condition then all overlying compartments are automatically assumed to
            # meet the condition. Thus in subsequent conditions we have to be careful
            # to ensure that True elements in 'cond4' do not also belong to 'cond5',
            # 'cond6' or 'cond7'. We use numpy.logical_not(...) for this purpose.
            cond5 = (self.var.th_fc_comp >= self.var.th_s_comp) & np.logical_not(cond4)
            cond6 = (zMid >= zGW_comp) & np.logical_not(cond4 | cond5)
            # cond6 = (zMid >= self.var.zGW) & np.logical_not(cond4 | cond5)
            cond7 = np.logical_not(cond4 | cond5 | cond6)
            dV = self.var.th_s_comp - self.var.th_fc_comp
            dFC = (dV / (Xmax ** 2)) * ((zMid - (zGW_comp - Xmax)) ** 2)
            # dFC = (dV / (Xmax ** 2)) * ((zMid - (self.var.zGW - Xmax)) ** 2)

            self.var.th_fc_adj[cond4] = self.var.th_fc_comp[cond4]
            self.var.th_fc_adj[cond5] = self.var.th_fc_comp[cond5]
            self.var.th_fc_adj[cond6] = self.var.th_s_comp[cond6]
            self.var.th_fc_adj[cond7] = self.var.th_fc_comp[cond7] + dFC[cond7]

        else:
            self.var.zGW = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon)) * -999
            self.var.WTinSoil = np.full((self.var.nCrop, self.var.nLat, self.var.nLon), False)
            self.var.th_fc_adj = np.copy(self.var.th_fc_comp)
