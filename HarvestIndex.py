#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from crop_growth_funs import *

import logging
logger = logging.getLogger(__name__)

class HarvestIndex(object):
    def __init__(self, HarvestIndex_variable):
        self.var = HarvestIndex_variable

    def initial(self):
        pass
    
    def dynamic(self):        
        """Function to simulate build up of harvest index"""
        
        # Get reference harvest index on current day
        HIi = np.copy(self.var.HIref)

        # Calculate harvest index
        # #######################
        
        HIadj = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        cond1 = (self.var.GrowingSeasonIndex & self.var.YieldForm & (self.var.HIt >= 0))

        # Root/tuber or fruit/grain crops
        cond11 = (cond1 & ((self.var.CropType == 2) | (self.var.CropType == 3)))

        # Determine adjustment for water stress before anthesis
        HI_adj_pre_anthesis(self.var.GrowingSeasonIndex, self.var.CropType, self.var.PreAdj, self.var.Fpre, self.var.YieldForm, self.var.HIt, self.var.B, self.var.B_NS, self.var.dHI_pre, self.var.CC)

        # Adjustment only for fruit/grain crops
        HImax = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))  # TODO: is this in the right place?
        cond112 = (cond11 & (self.var.CropType == 3))
        Ks = np.minimum(self.var.Ksw_Pol, self.var.Kst_PolC, self.var.Kst_PolH)
        HI_adj_pollination(self.var.GrowingSeasonIndex, self.var.CropType, self.var.Fpol, self.var.YieldForm, self.var.HIt, self.var.FloweringCD, self.var.CC, self.var.CCmin, Ks, self.var.exc)        
        HImax[cond112] = (self.var.Fpol * self.var.HI0)[cond112]
        cond113 = (cond11 & np.logical_not(cond112))
        HImax[cond113] = self.var.HI0[cond113]
                               
        # Determine adjustments for post-anthesis water stress
        harvest_index_adj_post_anthesis(self.var.GrowingSeasonIndex, self.var.CropType, self.var.Fpost, self.var.Fpre, self.var.sCor1, self.var.sCor2, self.var.fpost_upp, self.var.fpost_dwn, self.var.HIt, self.var.a_HI, self.var.b_HI, self.var.HIstartCD, self.var.HIendCD, self.var.YieldForm, self.var.YldFormCD, self.var.DAP, self.var.CC, self.var.CanopyDevEndCD, self.var.DelayedCDs, self.var.Ksw_Exp, self.var.Ksw_Sto)

        # Limit HI to maximum allowable increase due to pre- and post-anthesis
        # water stress combinations
        HImult = np.ones((self.var.nRotation, self.var.nLat, self.var.nLon))
        HImult[cond11] = (self.var.Fpre * self.var.Fpost)[cond11]
        cond115 = (cond11 & (HImult > (1 + (self.var.dHI0 / 100))))
        HImult[cond115] = (1 + (self.var.dHI0 / 100))[cond115]

        # Determine harvest index on current day, adjusted for stress effects
        cond116 = (cond11 & (HImax >= HIi))
        HIadj[cond116] = (HImult * HIi)[cond116]        
        cond117 = (cond11 & np.logical_not(cond116))
        HIadj[cond117] = (HImult * HImax)[cond117]

        # Leafy vegetable crops - no adjustment, harvest index equal to
        # reference value for current day
        cond12 = (cond1 & (self.var.CropType == 1))
        HIadj[cond12] = HIi[cond12]

        # Otherwise no build-up of harvest index if outside yield formation
        # period
        cond2 = (self.var.GrowingSeasonIndex & np.logical_not(cond1))
        HIi[cond2] = self.var.HI[cond2]
        HIadj[cond2] = self.var.HIadj[cond2]

        # Store final values for current time step
        self.var.HI[self.var.GrowingSeasonIndex] = HIi[self.var.GrowingSeasonIndex]
        self.var.HIadj[self.var.GrowingSeasonIndex] = HIadj[self.var.GrowingSeasonIndex]

        # No harvestable crop outside of a growing season
        self.var.HI[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.HIadj[np.logical_not(self.var.GrowingSeasonIndex)] = 0
