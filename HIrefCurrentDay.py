#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class HIrefCurrentDay(object):

    def __init__(self, HIrefCurrentDay_variable):
        self.var = HIrefCurrentDay_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.YieldForm = np.copy(arr_zeros)
        self.var.HIref = np.copy(arr_zeros)
        self.var.HIt = np.copy(arr_zeros)
        self.var.PctLagPhase = np.copy(arr_zeros)

    def reset_initial_conditions(self):
        self.var.PctLagPhase[self.var.GrowingSeasonDayOne] = 0
        
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

        # Yet to reach time for HI build-up
        cond1 = (self.var.GrowingSeasonIndex & (self.var.HIt <= 0))
        self.var.HIref[cond1] = 0
        self.var.PctLagPhase[cond1] = 0

        cond2 = (self.var.GrowingSeasonIndex & np.logical_not(cond1))

        # HI cannot develop further as canopy is too small (no need to do
        # anything here as self.var.HIref doesn't change)
        cond21 = (cond2 & (self.var.CCprev <= (self.var.CCmin * self.var.CCx)))
        cond22 = (cond2 & np.logical_not(cond21))

        # If crop type is leafy vegetable or root/tuber then proceed with
        # logistic growth (i.e. no linear switch)
        cond221 = (cond22 & ((self.var.CropType == 1) | (self.var.CropType == 2)))
        self.var.PctLagPhase[cond221] = 100
        HIref_divd = (self.var.HIini * self.var.HI0)
        HIref_divs = (self.var.HIini + (self.var.HI0 - self.var.HIini) * np.exp(-self.var.HIGC * self.var.HIt))
        HIref = np.divide(HIref_divd, HIref_divs, out=np.zeros_like(HIref_divs), where=HIref_divs!=0)
        self.var.HIref[cond221] = HIref[cond221]
        # print 'HIref1   %f' % self.var.HIref[0,0,0]
        
        # Harvest index approaching maximum limit
        cond2211 = (cond221 & (self.var.HIref >= (0.9799 * self.var.HI0)))
        self.var.HIref[cond2211] = self.var.HI0[cond2211]
        # print 'HIref2   %f' % self.var.HIref[0,0,0]

        cond222 = (cond22 & (self.var.CropType == 3))
        # Not yet reached linear switch point, therefore proceed with logistic
        # build-up
        cond2221 = (cond222 & (self.var.HIt < self.var.tLinSwitch))

        self.var.PctLagPhase[cond2221] = (100 * np.divide(self.var.HIt, self.var.tLinSwitch, out=np.zeros_like(self.var.HIt), where=self.var.tLinSwitch!=0))[cond2221]
        self.var.HIref[cond2221] = HIref[cond2221]
        # print 'HIref3   %f' % self.var.HIref[0,0,0]

        cond2222 = (cond222 & np.logical_not(cond2221))
        self.var.PctLagPhase[cond2222] = 100
        HIref_divd = (self.var.HIini * self.var.HI0)
        HIref_divs = (self.var.HIini + (self.var.HI0 - self.var.HIini) * np.exp(-self.var.HIGC * self.var.tLinSwitch))
        HIref = np.divide(HIref_divd, HIref_divs, out=np.zeros_like(HIref_divs), where=HIref_divs!=0)
        self.var.HIref[cond2222] = HIref[cond2222]

        # Calculate reference harvest index for current day (total = logistic + linear)
        self.var.HIref[cond2222] += (self.var.dHILinear * (self.var.HIt - self.var.tLinSwitch))[cond2222]
        # print 'HIref4   %f' % self.var.HIref[0,0,0]

        # print 'HIref      %f' % self.var.HIref[0,0,0]
        # print 'HIini      %f' % self.var.HIini[0,0,0]
        # print 'HI0        %f' % self.var.HI0[0,0,0]
        # print 'HIGC       %f' % self.var.HIGC[0,0,0]
        # print 'dHILinear  %f' % self.var.dHILinear[0,0,0]
        # print 'HIt        %f' % self.var.HIt[0,0,0]
        # print 'tLinSwitch %f' % self.var.tLinSwitch[0,0,0]
        
        # Limit self.var.HIref and round off computed value
        cond223 = (cond22 & (self.var.HIref > self.var.HI0))
        self.var.HIref[cond223] = self.var.HI0[cond223]
        cond224 = (cond22 & (self.var.HIref <= (self.var.HIini + 0.004)))
        self.var.HIref[cond224] = 0
        cond225 = (cond22 & ((self.var.HI0 - self.var.HIref) < 0.004))
        self.var.HIref[cond225] = self.var.HI0[cond225]

        # Reference harvest index is zero outside growing season
        self.var.HIref[np.logical_not(self.var.GrowingSeasonIndex)] = 0
