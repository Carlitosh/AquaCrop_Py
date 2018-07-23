#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class Evapotranspiration(object):
    """Class to represent combined soil evaporation and plant
    transpiration
    """
    def __init__(self, Evapotranspiration_variable):
        self.var = Evapotranspiration_variable

class AQEvapotranspiration(Evapotranspiration):
    
    def initial(self):
        self.var.ETpot = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

    def dynamic(self):
        self.var.ETpot = self.var.Epot + self.var.Tpot

class FAO56Evapotranspiration(Evapotranspiration):
    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.ETpot = np.copy(arr_zeros)
        self.var.ETpotCum = np.copy(arr_zeros)
        self.var.ETactCum = np.copy(arr_zeros)

    def dynamic(self):
        """Update Evapotranspiration for current day"""

        # Obtain crop coefficient for current growth stage
        L_day = np.stack((self.var.L_ini_day,
                          self.var.L_dev_day,
                          self.var.L_mid_day,
                          self.var.L_late_day), axis=0)        
        L_day = np.cumsum(L_day, axis=0)
        cond1 = self.var.GrowthStage == 1
        cond2 = self.var.GrowthStage == 2
        cond3 = self.var.GrowthStage == 3
        cond4 = self.var.GrowthStage == 4
        Kc = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        Kc[cond1] = self.var.Kc_ini[cond1]

        ini_to_mid_gradient = np.divide((self.var.Kc_mid - self.var.Kc_ini), self.var.L_dev_day, out=np.zeros_like(self.var.L_dev_day), where=self.var.L_dev_day!=0)

        Kc[cond2] = (self.var.Kc_ini + (ini_to_mid_gradient * (self.var.DAP - L_day[0,:])))[cond2]
        Kc[cond3] = self.var.Kc_mid[cond3]

        mid_to_end_gradient = np.divide((self.var.Kc_end - self.var.Kc_mid), self.var.L_late_day, out=np.zeros_like(self.var.L_late_day), where=self.var.L_late_day!=0)
        
        Kc[cond4] = (self.var.Kc_mid + (mid_to_end_gradient * (self.var.DAP - L_day[2,:])))[cond4]
        Kc[np.logical_not(self.var.GrowingSeasonIndex)] = 0.5  # Global Crop Water Model

        # Calculate maximum crop evapotranspiration
        self.var.ETpot = np.broadcast_to(self.var.referencePotET, (self.var.nCrop, self.var.nLat, self.var.nLon)) * Kc

        # Compute actual evapotranspiration
        p_std = np.copy(self.var.p_std)
        p_std[np.logical_not(self.var.GrowingSeasonIndex)] = 0.55  # Global Crop Water Model
        p = p_std + 0.04 * (5 - self.var.ETpot)  # Equation 31, root zone depletion factor
        Ks = (self.var.TAW - self.var.Dr) / ((1 - p) * self.var.TAW)  # Equation 30, water stress coefficient
        Ks = np.clip(Ks, 0, 1)
        self.var.ETact = self.var.ETpot * Ks  # Equation 29

        # Adjust water content to account for ET

        # Limit ETact to maximum available water in the root zone
        Wr_comp = self.var.RootFact * 1000 * self.var.th * self.var.dz_xy
        WrWP_comp = self.var.RootFact * 1000 * self.var.th_wp_comp * self.var.dz_xy        
        max_available_water = np.clip(Wr_comp - WrWP_comp, 0, None)
        self.var.ETact = np.clip(self.var.ETact, 0, np.sum(max_available_water, axis=1))
        
        # Extract water proportionally - TODO: this should be improved!!!
        f = (max_available_water / np.broadcast_to(np.sum(max_available_water, axis=1)[:,None,:,:], max_available_water.shape))
        ToExtract = np.broadcast_to(self.var.ETact[:,None,...], self.var.th.shape) * f
        ThToExtract = ((ToExtract / 1000) / self.var.dz_xy)
        self.var.th -= ThToExtract

        # Accumulate ETpot and ETact
        self.var.ETpotCum += self.var.ETpot
        self.var.ETactCum += self.var.ETact
