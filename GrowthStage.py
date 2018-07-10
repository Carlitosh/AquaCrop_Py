#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class GrowthStage(object):
    def __init__(self, GrowthStage_variable):
        self.var = GrowthStage_variable

    def initial(self):
        pass

    def dynamic(self):

        if self.var.CalendarType == 1:
            tAdj = self.var.DAP - self.var.DelayedCDs
        elif self.var.CalendarType == 2:
            tAdj = self.var.GDDcum - self.var.DelayedGDDs

        # Update growth stage
        cond1 = (self.var.GrowingSeasonIndex & (tAdj <= self.var.Canopy10Pct))
        cond2 = (self.var.GrowingSeasonIndex & np.logical_not(cond1) & (tAdj <= self.var.MaxCanopy))
        cond3 = (self.var.GrowingSeasonIndex & np.logical_not(cond1 | cond2) & (tAdj <= self.var.Senescence))
        cond4 = (self.var.GrowingSeasonIndex & np.logical_not(cond1 | cond2 | cond3) & (tAdj > self.var.Senescence))

        self.var.GrowthStage[cond1] = 1
        self.var.GrowthStage[cond2] = 2
        self.var.GrowthStage[cond3] = 3
        self.var.GrowthStage[cond4] = 4
        self.var.GrowthStage[np.logical_not(self.var.GrowingSeasonIndex)] = 0
