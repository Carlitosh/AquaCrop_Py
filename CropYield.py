#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class CropYield(object):
    def __init__(self, CropYield_variable):
        self.var = CropYield_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.var.Y = arr_zeros
        
    def dynamic(self):
        """Function to calculate crop yield"""
        cond1 = self.var.GrowingSeasonIndex
        self.var.Y[cond1] = ((self.var.B / 100) * self.var.HIadj)[cond1]
        cond11 = (cond1 & (((self.var.CalendarType == 1) & ((self.var.DAP - self.var.DelayedCDs) >= self.var.Maturity)) | ((self.var.CalendarType == 2) & ((self.var.GDDcum - self.var.DelayedGDDs) >= self.var.Maturity))))
        self.var.CropMature[cond11] = True
        self.var.Y[np.logical_not(self.var.GrowingSeasonIndex)] = 0
