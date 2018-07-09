#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from crop_growth_funs import *

import logging
logger = logging.getLogger(__name__)

class BiomassAccumulation(object):

    def __init__(self, BiomassAccumulation_variable):
        self.var = BiomassAccumulation_variable

    def initial(self):
        pass
    
    def dynamic(self):
        """Function to calculate biomass accumulation"""

        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        
        et0 = self.var.referencePotET[None,:,:] * np.ones((self.var.nRotation))[:,None,None]
        self.var.temperature_stress_module.dynamic()
        WPadj = adjust_WP_for_reproductive_stage(
            self.var.GrowingSeasonIndex, self.var.CropType, self.var.HIref, self.var.HIt, self.var.PctLagPhase, self.var.Determinant, self.var.YldFormCD, self.var.WP, self.var.WPy)
        
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
