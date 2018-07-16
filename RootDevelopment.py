#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class RootDevelopment(object):
    def __init__(self, RootDevelopment_variable):
        self.var = RootDevelopment_variable

class AQRootDevelopment(RootDevelopment):
    def initial(self):
        self.var.rCor = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.Zroot = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

    def reset_initial_conditions(self):
        self.var.rCor[self.var.GrowingSeasonDayOne] = 1
        self.var.Zroot[self.var.GrowingSeasonDayOne] = self.var.Zmin[self.var.GrowingSeasonDayOne]
        
    def potential_root_development(self, tAdj):
        dims = self.var.GrowingSeasonIndex.shape
        nr, nlat, nlon = dims[0], dims[1], dims[2]
        Zr = np.zeros((nr, nlat, nlon))
        Zini = self.var.Zmin * (self.var.PctZmin / 100)
        t0 = np.round(self.var.Emergence / 2)
        tmax = np.copy(self.var.MaxRooting)
        cond6 = (self.var.GrowingSeasonIndex & (tAdj >= tmax))
        Zr[cond6] = self.var.Zmax[cond6]
        cond7 = (self.var.GrowingSeasonIndex & np.logical_not(cond6) & (tAdj <= t0))
        Zr[cond7] = Zini[cond7]
        cond8 = (self.var.GrowingSeasonIndex & (np.logical_not(cond6 | cond7)))
        X_divd = (tAdj - t0).astype('float64')
        X_divs = (tmax - t0).astype('float64')
        X = np.divide(X_divd, X_divs, out=np.zeros_like(X_divs), where=X_divs!=0)
        Zr_exp = np.divide(1, self.var.fshape_r, out=np.zeros_like(self.var.fshape_r), where=cond8)
        Zr_pow = np.power(X, Zr_exp, out=np.zeros_like(Zr), where=cond8)
        Zr[cond8] = (Zini + (self.var.Zmax - Zini) * Zr_pow)[cond8]
        cond9 = (self.var.GrowingSeasonIndex & (Zr < self.var.Zmin))
        Zr[cond9] = self.var.Zmin[cond9]
        return Zr
        
    def dynamic(self):
        """Function to calculate root zone expansion"""

        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()
            
        dims = self.var.GrowingSeasonIndex.shape
        nr, nlat, nlon = dims[0], dims[1], dims[2]

        # Adjust time for any delayed development
        if self.var.CalendarType == 1:
            tAdj = (self.var.DAP - self.var.DelayedCDs)
        elif self.var.CalendarType == 2:
            tAdj = (self.var.GDDcum - self.var.DelayedGDDs)

        # Calculate root expansion
        if self.var.CalendarType == 1:
            tOld = (tAdj - 1)
        elif self.var.CalendarType == 2:
            tOld = (tAdj - self.var.GDD)

        tAdj[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        tOld[np.logical_not(self.var.GrowingSeasonIndex)] = 0

        # Potential root depth on previous day
        ZrOld = self.potential_root_development(tOld)
        
        # Potential root depth on current day
        Zr = self.potential_root_development(tAdj)

        # Determine rate of change, adjust for any stomatal water stress
        dZr = Zr - ZrOld
        cond10 = (self.var.GrowingSeasonIndex & (self.var.TrRatio < 0.9999))
        cond101 = (cond10 & (self.var.fshape_ex >= 0))        
        dZr[cond101] = (dZr * self.var.TrRatio)[cond101]
        cond102 = (cond10 & np.logical_not(cond101))
        fAdj_divd = (np.exp(self.var.TrRatio * self.var.fshape_ex) - 1)
        fAdj_divs = (np.exp(self.var.fshape_ex) - 1)
        fAdj = np.divide(fAdj_divd, fAdj_divs, out=np.zeros_like(Zr), where=fAdj_divs!=0)
        dZr[cond102] = (dZr * fAdj)[cond102]

        # Adjust root expansion for failure to germinate (roots cannot expand
        # if crop has not germinated)
        dZr[np.logical_not(self.var.Germination)] = 0

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

class FAO56RootDevelopment(RootDevelopment):
    def initial(self):
        self.var.Zroot = self.var.Zmin

    def dynamic(self):
        self.var.Zroot[self.var.GrowingSeasonIndex] = self.var.Zmin[self.var.GrowingSeasonIndex]
        self.var.Zroot[np.logical_not(self.var.GrowingSeasonIndex)] = 1.  # Global Crop Water Model
