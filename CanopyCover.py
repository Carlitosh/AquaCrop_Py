#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from crop_growth_funs import *

import logging
logger = logging.getLogger(__name__)

class CanopyCover(object):
    """Class to represent canopy growth/decline"""

    def __init__(self, CanopyCover_variable):
        self.var = CanopyCover_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.CCprev = arr_zeros
        
    def dynamic(self):    
        # Preallocate some variables
        CCxAdj = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        CDCadj = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        CCsen = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))

        # Store initial condition
        self.var.CCprev = np.copy(self.var.CC)
        self.var.CC_NSprev = np.copy(self.var.CC_NS)
        
        # Calculate root zone water content, and determine if water stress is
        # occurring
        # self.var.root_zone_water()
        # self.var.var.water_stress.dynamic(beta = True)

        # Get canopy cover growth over time
        if self.var.CalendarType == 1:
            tCC = np.copy(self.var.DAP)
            dtCC = np.ones((self.var.nRotation, self.var.nLat, self.var.nLon))
            tCCadj = self.var.DAP - self.var.DelayedCDs
        elif self.var.CalendarType == 2:
            tCC = np.copy(self.var.GDDcum)
            dtCC = np.copy(self.var.GDD)
            tCCadj = self.var.GDDcum - self.var.DelayedGDDs

        # Canopy development (potential)
        # ##############################
        self.var.CC_NS, self.var.CCxAct_NS, self.var.CCxW_NS = potential_canopy_development(
            self.var.GrowingSeasonIndex,
            self.var.Emergence, self.var.Maturity, self.var.Senescence,
            tCC, dtCC,
            self.var.CanopyDevEnd,
            self.var.CC0, self.var.CGC, self.var.CCx, self.var.CDC,
            self.var.CC_NS, self.var.CCxAct_NS, self.var.CCxW_NS)

        # Canopy development (actual)
        # ###########################
        self.var.CC, self.var.CC0adj, self.var.CCxAct = actual_canopy_development(
            self.var.GrowingSeasonIndex,
            self.var.Emergence, self.var.Maturity, self.var.Senescence,
            tCC, tCCadj, dtCC,
            self.var.CanopyDevEnd,
            self.var.CC0, self.var.CC0adj, self.var.CGC, self.var.CCx, self.var.CDC, self.var.Ksw_Exp, self.var.CC, self.var.CCprev, self.var.CCxAct)
        
        # Check for crop growth termination: if the following conditions are
        # met, the crop has died
        cond6 = (self.var.GrowingSeasonIndex & ((tCCadj <= self.var.Maturity) & (tCCadj >= self.var.Emergence) & (tCCadj > self.var.CanopyDevEnd)))
        cond63 = (cond6 & ((self.var.CC < 0.001) & np.logical_not(self.var.CropDead)))
        self.var.CC[cond63] = 0
        self.var.CropDead[cond63] = True

        # Canopy senescence due to water stress (actual)
        # ##############################################
        cond7 = (self.var.GrowingSeasonIndex & (tCCadj >= self.var.Emergence))

        # Check for early canopy senescence starting/continuing due to severe
        # water stress
        cond71 = (cond7 & ((tCCadj < self.var.Senescence) | (self.var.tEarlySen > 0)))

        # Early canopy senescence
        cond711 = (cond71 & (self.var.Ksw_Sen < 1))
        self.var.PrematSenes[cond711] = True

        # No prior early senescence
        cond7111 = (cond711 & (self.var.tEarlySen == 0))
        self.var.CCxEarlySen[cond7111] = self.var.CCprev[cond7111]

        # Increment early senescence GDD counter
        self.var.tEarlySen[cond711] += dtCC[cond711]

        update_CC_after_senescence(
            self.var.GrowingSeasonIndex, self.var.CC, self.var.CCprev, self.var.CC0, self.var.CC0adj, self.var.CGC, self.var.CDC, self.var.CCx, self.var.CCxAct, self.var.CCxW, self.var.CCxEarlySen, self.var.Ksw_Sen, self.var.tEarlySen, tCCadj, dtCC, self.var.Emergence, self.var.Senescence, self.var.PrematSenes, self.var.CropDead)

        # ##################################
        # adjust for micro-advective effects
        
        # Check to ensure potential CC is not slightly lower than actual
        cond8 = (self.var.GrowingSeasonIndex & (self.var.CC_NS < self.var.CC))
        self.var.CC_NS[cond8] = self.var.CC[cond8]

        cond81 = (cond8 & (tCC < self.var.CanopyDevEnd))
        self.var.CCxAct_NS[cond81] = self.var.CC_NS[cond81]

        # Actual (with water stress)
        self.var.CCadj[self.var.GrowingSeasonIndex] = ((1.72 * self.var.CC) - (self.var.CC ** 2) + (0.3 * (self.var.CC ** 3)))[self.var.GrowingSeasonIndex]

        # Potential (without water stress)
        self.var.CCadj_NS[self.var.GrowingSeasonIndex] = ((1.72 * self.var.CC_NS) - (self.var.CC_NS ** 2) + (0.3 * (self.var.CC_NS ** 3)))[self.var.GrowingSeasonIndex]

        # TODO: I don't think it is necessary to set these variables to zero here
        # No canopy outside growing season - set values to zero
        self.var.CC[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCadj[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CC_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCadj_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCxW[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCxAct[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCxW_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCxAct_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        
