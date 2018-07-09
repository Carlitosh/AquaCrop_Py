#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from crop_growth_funs import *

import logging
logger = logging.getLogger(__name__)

class HIrefCurrentDay(object):

    def __init__(self, HIrefCurrentDay_variable):
        self.var = HIrefCurrentDay_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.var.YieldForm = arr_zeros
        self.var.HIref = arr_zeros
        self.var.HIt = arr_zeros
        
    def dynamic(self):
        """Function to calculate reference (no adjustment for stress 
        effects) harvest index on current day
        """
        # Check if in yield formation period
        if self.var.CalendarType == 1:
            tAdj = self.var.DAP - self.var.DelayedCDs
        elif self.var.CalendarType == 2:
            tAdj = self.var.GDDcum - self.var.DelayedGDDs
        self.var.YieldForm = (self.var.GrowingSeasonIndex & (tAdj > self.var.HIstart))
        
        # Get time for harvest index calculation
        self.var.HIt = self.var.DAP - self.var.DelayedCDs - self.var.HIstartCD - 1

        harvest_index_ref_current_day(
            self.var.GrowingSeasonIndex, self.var.CropType,
            tAdj,
            self.var.HIt, self.var.HIref, self.var.HIini, self.var.HIGC, self.var.HI0, self.var.PctLagPhase,
            self.var.CCprev, self.var.CCmin, self.var.CCx, self.var.tLinSwitch, self.var.dHILinear)
