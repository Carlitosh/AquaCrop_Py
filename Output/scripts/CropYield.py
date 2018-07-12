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
        self.var.Y = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        
    def dynamic(self):
        """Function to calculate crop yield"""
        cond1 = self.var.GrowingSeasonIndex
        self.var.Y[cond1] = ((self.var.B / 100) * self.var.HIadj)[cond1]
        cond11 = (cond1 & (((self.var.CalendarType == 1) & ((self.var.DAP - self.var.DelayedCDs) >= self.var.Maturity)) | ((self.var.CalendarType == 2) & ((self.var.GDDcum - self.var.DelayedGDDs) >= self.var.Maturity))))
        self.var.CropMature[cond11] = True
        self.var.Y[np.logical_not(self.var.GrowingSeasonIndex)] = 0

class FAO56CropYield(CropYield):
    
    def initial(self):
        self.var.Y = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))

    def dynamic(self):

        # TODO: only read data if the year has changed
        
        date = '%04i-%02i-%02i' %(self.var._modelTime.year, 1, 1)
        self.var.Yx = vos.netcdf2PCRobjClone(self.var.PotYieldFileNC,
                                             self.var.PotYieldVarName,
                                             date,
                                             useDoy = None,
                                             cloneMapFileName = self.var.cloneMap,
                                             LatitudeLongitude = True)
        cond1 = self.var._modelTime.doy == self.var.HarvestDateAdj
        self.var.Y[cond1] = (self.var.Yx * (1 - self.var.Ky * (1 - self.var.ETactCum / self.var.ETpotCum)))[cond1]
        self.var.Y[np.logical_not(cond1)] = 0
        print self.var.Y[0,0,0]
        
        
        
