#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from crop_growth_funs import *

import logging
logger = logging.getLogger(__name__)

class TemperatureStress(object):
    def __init__(self, TemperatureStress_variable):
        self.var = TemperatureStress_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.Kst_Bio = np.copy(arr_zeros)
        self.Kst_PolH = np.copy(arr_zeros)
        self.Kst_PolC = np.copy(arr_zeros)
        
    def dynamic(self):
        """Function to calculate temperature stress coefficients"""

        # Add rotation dimension to meteo vars
        tmin = self.var.tmin[None,:,:] * np.ones((self.var.nRotation))[:,None,None]
        tmax = self.var.tmax[None,:,:] * np.ones((self.var.nRotation))[:,None,None]
        
        # Calculate temperature stress coefficients affecting crop pollination
        self.var.Kst_Bio = temperature_stress_biomass(self.var.BioTempStress, self.var.GDD, self.var.GDD_up, self.var.GDD_lo)

        KsPol_up = 1
        KsPol_lo = 0.001
        self.var.Kst_PolH = temperature_stress_heat(tmax, self.var.Tmax_lo, self.var.Tmax_up, self.var.PolHeatStress, self.var.fshape_b, KsPol_up, KsPol_lo)
        self.var.Kst_PolC = temperature_stress_cold(tmin, self.var.Tmin_lo, self.var.Tmin_up, self.var.PolColdStress, self.var.fshape_b, KsPol_up, KsPol_lo)
