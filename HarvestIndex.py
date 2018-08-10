#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class HarvestIndex(object):
    def __init__(self, HarvestIndex_variable):
        self.var = HarvestIndex_variable

    def initial(self):
        arr_ones = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.Fpre = np.copy(arr_ones)
        self.var.Fpost = np.copy(arr_ones)
        self.var.fpost_dwn = np.copy(arr_ones)
        self.var.fpost_upp = np.copy(arr_ones)
        self.var.Fpol = np.copy(arr_zeros)
        self.var.sCor1 = np.copy(arr_zeros)
        self.var.sCor2 = np.copy(arr_zeros)
        self.var.HI = np.copy(arr_zeros)
        self.var.HIadj = np.copy(arr_zeros)
        self.var.PreAdj = np.copy(arr_zeros.astype(bool))
        
    def reset_initial_conditions(self):
        self.var.Fpre[self.var.GrowingSeasonDayOne] = 1
        self.var.Fpost[self.var.GrowingSeasonDayOne] = 1
        self.var.fpost_dwn[self.var.GrowingSeasonDayOne] = 1
        self.var.fpost_upp[self.var.GrowingSeasonDayOne] = 1
        self.var.Fpol[self.var.GrowingSeasonDayOne] = 0
        self.var.sCor1[self.var.GrowingSeasonDayOne] = 0
        self.var.sCor2[self.var.GrowingSeasonDayOne] = 0
        self.var.HI[self.var.GrowingSeasonDayOne] = 0
        self.var.HIadj[self.var.GrowingSeasonDayOne] = 0
        self.var.PreAdj[self.var.GrowingSeasonDayOne] = False

    def HI_adj_pollination(self):
        """Function to calculate adjustment to harvest index for 
        failure of pollination due to water or temperature stress
        """
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        FracFlow = np.copy(arr_zeros)
        t1 = np.copy(arr_zeros)
        t2 = np.copy(arr_zeros)
        F1 = np.copy(arr_zeros)
        F2 = np.copy(arr_zeros)
        F = np.copy(arr_zeros)

        cond0 = (self.var.GrowingSeasonIndex & self.var.YieldForm & (self.var.CropType == 3) & (self.var.HIt > 0) & (self.var.HIt <= self.var.FloweringCD))

        # Fractional flowering on previous day
        # cond1 = (HIt > 0)
        t1[cond0] = self.var.HIt[cond0] - 1
        cond11 = (cond0 & (t1 > 0))
        t1pct = 100 * np.divide(t1, self.var.FloweringCD, out=np.copy(arr_zeros), where=self.var.FloweringCD!=0)
        t1pct = np.clip(t1pct, 0, 100)
        F1[cond11] = (0.00558 * np.exp(0.63 * np.log(t1pct, out=np.copy(arr_zeros), where=t1pct>0)) - (0.000969 * t1pct) - 0.00383)[cond11]
        F1 = np.clip(F1, 0, None)

        # Fractional flowering on current day
        t2[cond0] = self.var.HIt[cond0]
        cond12 = (cond0 & (t2 > 0))
        t2pct = 100 * np.divide(t2, self.var.FloweringCD, out=np.copy(arr_zeros), where=self.var.FloweringCD!=0)
        t2pct = np.clip(t2pct, 0, 100)
        F2[cond12] = (0.00558 * np.exp(0.63 * np.log(t2pct, out=np.copy(arr_zeros), where=t2pct>0)) - (0.000969 * t2pct) - 0.00383)[cond12]
        F2 = np.clip(F2, 0, None)

        # Weight values
        cond13 = (cond0 & (np.abs(F1 - F2) >= 0.0000001))
        F[cond13] = (100 * np.divide(((F1 + F2) / 2), self.var.FloweringCD, out=np.copy(arr_zeros), where=self.var.FloweringCD!=0))[cond13]
        FracFlow[cond13] = F[cond13]

        # Calculate pollination adjustment for current day
        dFpol = np.copy(arr_zeros)##np.zeros((self.nCrop, self.nLat, self.nLon))
        cond2 = (cond0 & (self.var.CC >= self.var.CCmin))
        Ks = np.minimum(self.var.Ksw_Pol, self.var.Kst_PolC, self.var.Kst_PolH)
        dFpol[cond2] = (Ks * FracFlow * (1 + (self.var.exc / 100)))[cond2]

        # Calculate pollination adjustment to dateppp
        self.var.Fpol += dFpol
        self.var.Fpol = np.clip(self.var.Fpol, None, 1)

    def HI_adj_pre_anthesis(self):
        """Function to calculate adjustment to harvest index for 
        pre-anthesis water stress
        """
        cond0 = (self.var.GrowingSeasonIndex & self.var.YieldForm & (self.var.HIt >= 0) & ((self.var.CropType == 2) | (self.var.CropType == 3)) & np.logical_not(self.var.PreAdj))
        self.var.PreAdj[cond0] = True

        # Calculate adjustment
        Br = np.divide(self.var.B, self.var.B_NS, out=np.zeros_like(self.var.B_NS), where=self.var.B_NS!=0)
        Br_range = np.log(self.var.dHI_pre, out=np.zeros_like(self.var.dHI_pre), where=self.var.dHI_pre>0) / 5.62
        Br_upp = 1
        Br_low = 1 - Br_range
        Br_top = Br_upp - (Br_range / 3)

        # Get biomass ratio
        ratio_low_divd = (Br - Br_low)
        ratio_low_divs = (Br_top - Br_low)
        ratio_low = np.divide(ratio_low_divd, ratio_low_divs, out=np.zeros_like(ratio_low_divs), where=ratio_low_divs!=0)
        ratio_upp_divd = (Br - Br_top)
        ratio_upp_divs = (Br_upp - Br_top)
        ratio_upp = np.divide(ratio_upp_divd, ratio_upp_divs, out=np.zeros_like(ratio_upp_divs), where=ratio_upp_divs!=0)

        # Calculate adjustment factor
        cond1 = (cond0 & ((Br >= Br_low) & (Br < Br_top)))
        self.var.Fpre[cond1] = (1 + (((1 + np.sin((1.5 - ratio_low) * np.pi)) / 2) * (self.var.dHI_pre / 100)))[cond1]
        cond2 = (cond0 & np.logical_not(cond1) & ((Br > Br_top) & (Br <= Br_upp)))
        self.var.Fpre[cond2] = (1 + (((1 + np.sin((0.5 + ratio_upp) * np.pi)) / 2) * (self.var.dHI_pre / 100)))[cond2]
        cond3 = (cond0 & np.logical_not(cond1 | cond2))
        self.var.Fpre[cond3] = 1

        # No green canopy left at start of flowering so no harvestable crop
        # will develop
        cond3 = (cond0 & (self.var.CC <= 0.01))
        self.var.Fpre[cond3] = 0

    def HI_adj_post_anthesis(self):
        """Function to calculate adjustment to harvest index for 
        post-anthesis water stress
        """
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        cond0 = (self.var.GrowingSeasonIndex & self.var.YieldForm & (self.var.HIt > 0) & ((self.var.CropType == 2) | (self.var.CropType == 3)))

        # 1 Adjustment for leaf expansion
        tmax1 = self.var.CanopyDevEndCD - self.var.HIstartCD
        self.var.DAP -= self.var.DelayedCDs
        cond1 = (cond0 & (self.var.DAP <= (self.var.CanopyDevEndCD + 1)) & (tmax1 > 0) & (self.var.Fpre > 0.99) & (self.var.CC > 0.001) & (self.var.a_HI > 0))
        dCor = (1 + np.divide((1 - self.var.Ksw_Exp), self.var.a_HI, out=np.copy(arr_zeros), where=self.var.a_HI!=0))
        self.var.sCor1[cond1] += np.divide(dCor, tmax1, out=np.copy(arr_zeros), where=tmax1!=0)[cond1]
        DayCor = (self.var.DAP - 1 - self.var.HIstartCD)
        self.var.fpost_upp[cond1] = (np.divide(tmax1, DayCor, out=np.copy(arr_zeros), where=DayCor!=0) * self.var.sCor1)[cond1]

        # 2 Adjustment for stomatal closure
        tmax2 = np.copy(self.var.YldFormCD)
        cond2 = (cond0 & (self.var.DAP <= (self.var.HIendCD + 1)) & (tmax2 > 0) & (self.var.Fpre > 0.99) & (self.var.CC > 0.001) & (self.var.b_HI > 0))
        dCor = ((np.exp(0.1 * np.log(self.var.Ksw_Sto, out=np.zeros_like(self.var.Ksw_Sto), where=self.var.Ksw_Sto!=0))) * (1 - np.divide((1 - self.var.Ksw_Sto), self.var.b_HI, out=np.copy(arr_zeros), where=self.var.b_HI!=0)))
        self.var.sCor2[cond2] += np.divide(dCor, tmax2, out=np.copy(arr_zeros), where=tmax2!=0)[cond2]
        DayCor = (self.var.DAP - 1 - self.var.HIstartCD)
        self.var.fpost_dwn[cond2] = (np.divide(tmax2, DayCor, out=np.copy(arr_zeros), where=DayCor!=0) * self.var.sCor2)[cond2]

        # Determine total multiplier
        cond3 = (cond0 & (tmax1 == 0) & (tmax2 == 0))
        self.var.Fpost[cond3] = 1
        cond4 = (cond0 & np.logical_not(cond3))
        cond41 = (cond4 & (tmax2 == 0))
        self.var.Fpost[cond41] = self.var.fpost_upp[cond41]
        cond42 = (cond4 & (tmax1 <= tmax2) & np.logical_not(cond41))
        self.var.Fpost[cond42] = (self.var.fpost_dwn * np.divide(((tmax1 * self.var.fpost_upp) + (tmax2 - tmax1)), tmax2, out=np.copy(arr_zeros), where=tmax2!=0))[cond42]
        cond43 = (cond4 & np.logical_not(cond41 | cond42))
        self.var.Fpost[cond43] = (self.var.fpost_upp * np.divide(((tmax2 * self.var.fpost_dwn) + (tmax1 - tmax2)), tmax2, out=np.copy(arr_zeros), where=tmax2!=0))[cond43]
        
    def dynamic(self):        
        """Function to simulate build up of harvest index"""

        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()
            
        # Get reference harvest index on current day
        HIi = np.copy(self.var.HIref)

        # Calculate harvest index
        # #######################
        
        HIadj = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        cond1 = (self.var.GrowingSeasonIndex & self.var.YieldForm & (self.var.HIt >= 0))

        # Root/tuber or fruit/grain crops
        cond11 = (cond1 & ((self.var.CropType == 2) | (self.var.CropType == 3)))

        # Determine adjustment for water stress before anthesis
        self.HI_adj_pre_anthesis()
        
        # Adjustment only for fruit/grain crops
        HImax = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))  # TODO: is this in the right place?
        cond112 = (cond11 & (self.var.CropType == 3))
        Ks = np.minimum(self.var.Ksw_Pol, self.var.Kst_PolC, self.var.Kst_PolH)
        self.HI_adj_pollination()
        HImax[cond112] = (self.var.Fpol * self.var.HI0)[cond112]
        cond113 = (cond11 & np.logical_not(cond112))
        HImax[cond113] = self.var.HI0[cond113]
                               
        # Determine adjustments for post-anthesis water stress
        self.HI_adj_post_anthesis()

        # Limit HI to maximum allowable increase due to pre- and post-anthesis
        # water stress combinations
        HImult = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
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
