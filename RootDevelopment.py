#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from crop_growth_funs import *

import logging
logger = logging.getLogger(__name__)

class RootDevelopment(object):
    def __init__(self, RootDevelopment_variable):
        self.var = RootDevelopment_variable

    def initial(self):
        pass
    
    def dynamic(self):
        """Function to calculate root zone expansion"""

        dZr, Zr = root_development(
            self.var.GrowingSeasonIndex,
            self.var.CalendarType,
            self.var.DAP, self.var.DelayedCDs, self.var.DelayedGDDs, self.var.GDD, self.var.GDDcum, self.var.Zmin, self.var.Zmax, self.var.PctZmin, self.var.Emergence, self.var.MaxRooting, self.var.fshape_r, self.var.fshape_ex, self.var.TrRatio, self.var.Germination)
        
        # Get new rooting depth
        self.var.Zroot[self.var.GrowingSeasonIndex] = (self.var.Zroot + dZr)[self.var.GrowingSeasonIndex]

        # Adjust root depth if restrictive soil layer is present that limits
        # depth of root expansion
        cond11 = (self.var.GrowingSeasonIndex & (self.var.zRes > 0))
        cond111 = (cond11 & (self.var.Zroot > self.var.zRes))
        self.var.rCor[cond111] = np.divide(
            (2 * (self.var.Zroot / self.var.zRes) * ((self.var.SxTop + self.var.SxBot) / 2) - self.var.SxTop),
            self.var.SxBot, out=np.zeros_like(Zr), where=cond111)[cond111]        
        self.var.Zroot[cond111] = self.var.zRes[cond111]

        # Limit rooting depth if groundwater table is present (roots cannot
        # develop below the water table)
        if self.var.WaterTable:
            zGW = np.copy(self.var.zGW)
            cond12 = ((zGW > 0) & (self.var.Zroot > zGW))
            self.var.Zroot[cond12] = np.clip(zGW, self.var.Zmin, None)[cond12]

        # No root system outside of growing season
        self.var.Zroot[np.logical_not(self.var.GrowingSeasonIndex)] = 0
