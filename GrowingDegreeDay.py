#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class GrowingDegreeDay(object):
    """Class to represent growing degree days"""

    def __init__(self, GrowingDegreeDay_variable):
        self.var = GrowingDegreeDay_variable

    def initial(self):
        self.var.GDDcum = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.GDD = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

    def reset_initial_conditions(self):
        self.var.GDDcum[self.var.GrowingSeasonDayOne] = 0

    def growing_degree_day(self):
        tmax = self.var.tmax[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
        tmin = self.var.tmin[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
        # tmax = self.var.tmax[None,:,:] * np.ones((self.var.nRotation))[:,None,None]
        # tmin = self.var.tmin[None,:,:] * np.ones((self.var.nRotation))[:,None,None]
        if self.var.GDDmethod == 1:
            Tmean = ((tmax + tmin) / 2)
            Tmean = np.clip(Tmean, self.var.Tbase, self.var.Tupp)
        elif self.var.GDDmethod == 2:
            tmax = np.clip(tmax, self.var.Tbase, self.var.Tupp)
            tmin = np.clip(tmin, self.var.Tbase, self.var.Tupp)
            Tmean = ((tmax + tmin) / 2)
        elif self.var.GDDmethod == 3:
            tmax = np.clip(tmax, self.var.Tbase, self.var.Tupp)
            tmin = np.clip(tmin, None, self.var.Tupp)
            Tmean = ((tmax + tmin) / 2)
            Tmean = np.clip(Tmean, self.var.Tbase, None)
        self.var.GDD = (Tmean - self.var.Tbase)
        
    def dynamic(self):
        self.growing_degree_day()
        self.var.GDDcum[self.var.GrowingSeasonIndex] += self.var.GDD[self.var.GrowingSeasonIndex]
        self.var.GDDcum[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        
        
        

