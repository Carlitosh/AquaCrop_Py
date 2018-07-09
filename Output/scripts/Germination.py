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
