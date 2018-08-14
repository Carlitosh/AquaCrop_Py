#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np
import VirtualOS as vos

import logging
logger = logging.getLogger(__name__)

class CropYield(object):
    def __init__(self, CropYield_variable):
        self.var = CropYield_variable

class AQCropYield(CropYield):
    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.CropMature = np.copy(arr_zeros).astype(bool)
        self.var.Y = np.copy(arr_zeros)

    def reset_initial_conditions(self):
        self.var.CropMature[self.var.GrowingSeasonDayOne] = False
        
    def dynamic(self):
        """Function to calculate crop yield"""
        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()
            
        cond1 = self.var.GrowingSeasonIndex
        self.var.Y[cond1] = ((self.var.B / 100) * self.var.HIadj)[cond1]
        cond11 = (cond1 & (((self.var.CalendarType == 1) & ((self.var.DAP - self.var.DelayedCDs) >= self.var.Maturity)) | ((self.var.CalendarType == 2) & ((self.var.GDDcum - self.var.DelayedGDDs) >= self.var.Maturity))))
        self.var.CropMature[cond11] = True
        self.var.Y[np.logical_not(self.var.GrowingSeasonIndex)] = 0

class FAO56CropYield(CropYield):
    
    def initial(self):
        self.var.Y = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

    def dynamic(self):
        
        cond1 = self.var._modelTime.doy == self.var.HarvestDateAdj
        self.var.Y[cond1] = (self.var.Yx * (1 - self.var.Ky * (1 - self.var.ETactCum / self.var.ETpotCum)))[cond1]
        self.var.Y[np.logical_not(cond1)] = 0
        # print self.var._modelTime.currTime.timetuple().tm_yday
        # print self.var.GrowingSeasonIndex[:,10,10]
        # print np.max(self.var.Y, axis=(1,2))
        
