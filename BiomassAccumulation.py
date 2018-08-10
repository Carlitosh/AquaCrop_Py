#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class BiomassAccumulation(object):

    def __init__(self, BiomassAccumulation_variable):
        self.var = BiomassAccumulation_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.B = np.copy(arr_zeros)
        self.var.B_NS = np.copy(arr_zeros)

    def reset_initial_conditions(self):
        self.var.B[self.var.GrowingSeasonDayOne] = 0
        self.var.B_NS[self.var.GrowingSeasonDayOne] = 0
        
    def dynamic(self):
        """Function to calculate biomass accumulation"""

        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()
            
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        
        et0 = self.var.referencePotET[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
        self.var.temperature_stress_module.dynamic()

        fswitch = np.copy(arr_zeros)
        WPadj = np.copy(arr_zeros)
        cond1 = (self.var.GrowingSeasonIndex & (((self.var.CropType == 2) | (self.var.CropType == 3)) & (self.var.HIref > 0)))

        # Adjust WP for reproductive stage
        cond11 = (cond1 & (self.var.Determinant == 1))
        fswitch[cond11] = (self.var.PctLagPhase / 100)[cond11]
        cond12 = (cond1 & np.logical_not(cond11))
        cond121 = (cond12 < (self.var.YldFormCD / 3))
        fswitch[cond121] = np.divide(self.var.HIt.astype(np.float64), (self.var.YldFormCD.astype(np.float64) / 3.), out=np.zeros_like(self.var.YldFormCD.astype(np.float64)), where=self.var.YldFormCD!=0)[cond121]
        cond122 = (cond12 & np.logical_not(cond121))
        fswitch[cond122] = 1
        WPadj[cond1] = (self.var.WP * (1 - (1 - self.var.WPy / 100) * fswitch))[cond1]
        cond2 = (self.var.GrowingSeasonIndex & np.logical_not(cond1))
        WPadj[cond2] = self.var.WP[cond2]
        
        # Adjust WP for CO2 effects)
        WPadj *= self.var.fCO2 

        # Calculate biomass accumulation on current day
        dB_NS = WPadj * (self.var.TrPot_NS / et0) * self.var.Kst_Bio  # TODO: check correct TrPot is being used
        dB = WPadj * (self.var.TrAct / et0) * self.var.Kst_Bio  # TODO: check correct TrAct is being used

        # Update biomass accumulation
        self.var.B += dB
        self.var.B_NS += dB_NS

        # No biomass accumulation outside growing season
        self.var.B[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.B_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
