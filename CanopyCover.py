#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class CanopyCover(object):
    """Class to represent canopy growth/decline"""

    def __init__(self, CanopyCover_variable):
        self.var = CanopyCover_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.tEarlySen = np.copy(arr_zeros)
        self.var.CC = np.copy(arr_zeros)
        self.var.CCadj = np.copy(arr_zeros)
        self.var.CC_NS = np.copy(arr_zeros)
        self.var.CCadj_NS = np.copy(arr_zeros)
        self.var.CCxAct = np.copy(arr_zeros)
        self.var.CCxAct_NS = np.copy(arr_zeros)
        self.var.CCxW = np.copy(arr_zeros)
        self.var.CCxW_NS = np.copy(arr_zeros)
        self.var.CCxEarlySen = np.copy(arr_zeros)
        self.var.CCprev = np.copy(arr_zeros)
        self.var.PrematSenes = np.copy(arr_zeros.astype(bool))        
        self.var.CropDead = np.copy(arr_zeros.astype(bool))        
        self.var.CC0adj = np.copy(arr_zeros)

    def reset_initial_conditions(self):
        self.var.tEarlySen[self.var.GrowingSeasonDayOne] = 0
        self.var.CC[self.var.GrowingSeasonDayOne] = 0
        self.var.CCadj[self.var.GrowingSeasonDayOne] = 0
        self.var.CC_NS[self.var.GrowingSeasonDayOne] = 0
        self.var.CCadj_NS[self.var.GrowingSeasonDayOne] = 0
        self.var.CCxAct[self.var.GrowingSeasonDayOne] = 0
        self.var.CCxAct_NS[self.var.GrowingSeasonDayOne] = 0
        self.var.CCxW[self.var.GrowingSeasonDayOne] = 0
        self.var.CCxW_NS[self.var.GrowingSeasonDayOne] = 0
        self.var.CCxEarlySen[self.var.GrowingSeasonDayOne] = 0
        self.var.CCprev[self.var.GrowingSeasonDayOne] = 0
        self.var.PrematSenes[self.var.GrowingSeasonDayOne] = False
        self.var.CropDead[self.var.GrowingSeasonDayOne] = False
        self.var.CC0adj[self.var.GrowingSeasonDayOne] = self.var.CC0[self.var.GrowingSeasonDayOne]
        
    def canopy_cover_development(self, CC0, CCx, CGC, CDC, dt, Mode):
        """Function to calculate canopy cover development by end of the 
        current simulation day
        """
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        if Mode == 'Growth':
            CC = (CC0 * np.exp(CGC * dt))
            cond1 = (CC > (CCx / 2.))
            CC[cond1] = (CCx - 0.25 * np.divide(CCx, CC0, out=np.copy(arr_zeros), where=CC0!=0) * CCx * np.exp(-CGC * dt))[cond1]
            CC = np.clip(CC, None, CCx)
        elif Mode == 'Decline':
            CC = np.copy(arr_zeros)
            cond2 = (CCx >= 0.001)
            CC[cond2] = (CCx * (1. - 0.05 * (np.exp(dt * np.divide(CDC, CCx, out=np.copy(arr_zeros), where=CCx!=0)) - 1.)))[cond2]

        CC = np.clip(CC, 0, 1)
        return CC

    def canopy_cover_required_time(self, CC0, CCx, CGC, CDC, dt, tSum, Mode):
        """Function to find required time to reach CC at end of previous 
        day, given current CGC or CDC
        """
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        if Mode == 'CGC':
            CGCx = np.copy(arr_zeros)
            cond1 = (self.var.CCprev <= (CCx / 2))
            x = np.divide(self.var.CCprev, CC0, out=np.copy(arr_zeros), where=CC0!=0)
            CGCx_divd = np.log(x, out=np.copy(arr_zeros), where=x>0)
            CGCx_divs = tSum - dt
            CGCx[cond1] = np.divide(CGCx_divd, CGCx_divs, out=np.copy(arr_zeros), where=CGCx_divs!=0)[cond1]
            cond2 = np.logical_not(cond1)

            x1 = np.divide(0.25 * CCx * CCx, CC0, out=np.copy(arr_zeros), where=CC0!=0)
            x2 = CCx - self.var.CCprev
            x3 = np.divide(x1, x2, out=np.copy(arr_zeros), where=x2!=0)
            CGCx_divd = np.log(x3, out=np.copy(arr_zeros), where=x3>0)
            CGCx_divs = tSum - dt
            CGCx[cond2] = np.divide(CGCx_divd, CGCx_divs, out=np.copy(arr_zeros), where=CGCx_divs!=0)[cond2]
            tReq = (tSum - dt) * np.divide(CGCx, CGC, out=np.copy(arr_zeros), where=CGC!=0)
        elif Mode == 'CDC':
            
            x1 = np.divide(self.var.CCprev, CCx, out=np.copy(arr_zeros), where=CCx!=0)
            x2 = 1 + (1 - x1) / 0.05
            tReq_divd = np.log(x2, out=np.copy(arr_zeros), where=x2!=0)
            tReq_divs = np.divide(CDC, CCx, out=np.copy(arr_zeros), where=CCx!=0)
            tReq = np.divide(tReq_divd, tReq_divs, out=np.copy(arr_zeros), where=tReq_divs!=0)
            
        return tReq

    def adjust_CCx(self, CC0, CCx, CGC, CDC, dt, tSum, CanopyDevEnd):
        """Function to adjust CCx value for changes in CGC due to water 
        stress during the growing season
        """
        # Get time required to reach CC on previous day, then calculate
        # adjusted CCx
        tCCtmp = self.canopy_cover_required_time(CC0, CCx, CGC, CDC, dt, tSum, 'CGC')
        cond1 = (tCCtmp > 0)
        tCCtmp[cond1] += ((CanopyDevEnd - tSum) + dt)[cond1]
        CCxAdj = self.canopy_cover_development(CC0, CCx, CGC, CDC, tCCtmp, 'Growth')
        CCxAdj[np.logical_not(cond1)] = 0
        return CCxAdj
    
    def potential_canopy_development(self, tCC, dtCC):

        CC_NSprev = np.copy(self.var.CC_NS)

        # No canopy development before emergence/germination or after maturity
        cond1 = (self.var.GrowingSeasonIndex & ((tCC < self.var.Emergence) | (np.round(tCC) > self.var.Maturity)))
        self.var.CC_NS[cond1] = 0

        # Canopy growth can occur
        cond2 = (self.var.GrowingSeasonIndex & np.logical_not(cond1) & (tCC < self.var.CanopyDevEnd))

        # Very small initial CC as it is first day or due to senescence. In this
        # case assume no leaf expansion stress
        cond21 = (cond2 & (CC_NSprev <= self.var.CC0))
        self.var.CC_NS[cond21] = (self.var.CC0 * np.exp(self.var.CGC * dtCC))[cond21]

        # Canopy growing
        cond22 = (cond2 & np.logical_not(cond21))
        tmp_tCC = tCC - self.var.Emergence
        tmp_CC_NS = self.canopy_cover_development(self.var.CC0, self.var.CCx, self.var.CGC, self.var.CDC, tmp_tCC, 'Growth')
        self.var.CC_NS[cond22] = tmp_CC_NS[cond22]

        # Update maximum canopy cover size in growing season
        self.var.CCxAct_NS[cond2] = self.var.CC_NS[cond2]

        # No more canopy growth is possible or canopy in decline
        cond3 = (self.var.GrowingSeasonIndex & np.logical_not(cond1 | cond2) & (tCC > self.var.CanopyDevEnd))
        # Set CCx for calculation of withered canopy effects
        self.var.CCxW_NS[cond3] = self.var.CCxAct_NS[cond3]

        # Mid-season stage - no canopy growth, so do not update CC_NS
        cond31 = (cond3 & (tCC < self.var.Senescence))
        self.var.CC_NS[cond31] = CC_NSprev[cond31]
        self.var.CCxAct_NS[cond31] = self.var.CC_NS[cond31]

        # Late-season stage - canopy decline
        cond32 = (cond3 & np.logical_not(cond31))
        tmp_tCC = tCC - self.var.Emergence
        tmp_CC_NS = self.canopy_cover_development(self.var.CC0, self.var.CCx, self.var.CGC, self.var.CDC, tmp_tCC, 'Decline')
        self.var.CC_NS[cond32] = tmp_CC_NS[cond32]
        # return CC_NS, CCxAct_NS, CCxW_NS

    def actual_canopy_development(self, tCC, tCCadj, dtCC):

        # CCprev = np.copy(CC)

        # No canopy development before emergence/germination or after maturity
        cond4 = (self.var.GrowingSeasonIndex & ((tCCadj < self.var.Emergence) | (np.round(tCCadj) > self.var.Maturity)))
        self.var.CC[cond4] = 0

        # Otherwise, canopy growth can occur
        cond5 = (self.var.GrowingSeasonIndex & np.logical_not(cond4) & (tCCadj < self.var.CanopyDevEnd))
        cond51 = (cond5 & (self.var.CCprev <= self.var.CC0adj))

        # Very small initial CC as it is first day or due to senescence. In
        # this case, assume no leaf expansion stress
        self.var.CC[cond51] = (self.var.CC0adj * np.exp(self.var.CGC * dtCC))[cond51]

        # Canopy growing
        cond52 = (cond5 & np.logical_not(cond51))

        # Canopy approaching maximum size
        cond521 = (cond52 & (self.var.CCprev >= (0.9799 * self.var.CCx)))
        tmp_tCC = tCC - self.var.Emergence
        tmp_CC = self.canopy_cover_development(self.var.CC0, self.var.CCx, self.var.CGC, self.var.CDC, tmp_tCC, 'Growth')
        self.var.CC[cond521] = tmp_CC[cond521]
        self.var.CC0adj[cond521] = self.var.CC0[cond521]

        # Adjust canopy growth coefficient for leaf expansion water stress
        # effects
        cond522 = (cond52 & np.logical_not(cond521))
        CGCadj = self.var.CGC * self.var.Ksw_Exp

        # Adjust CCx for change in CGC
        cond5221 = (cond522 & (CGCadj > 0))
        CCxAdj = self.adjust_CCx(self.var.CC0adj, self.var.CCx, CGCadj, self.var.CDC, dtCC, tCCadj, self.var.CanopyDevEnd)
        # CCxAdj = self.adjust_CCx(self.var.CCprev, self.var.CC0adj, self.var.CCx, CGCadj, self.var.CDC, dtCC, tCCadj, self.var.CanopyDevEnd)
        cond52211 = (cond5221 & (CCxAdj > 0))
        # Approaching maximum canopy size
        cond522111 = (cond52211 & (np.abs(self.var.CCprev - self.var.CCx) < 0.00001))
        tmp_tCC = tCC - self.var.Emergence
        # tmp_CC = self.canopy_cover_development(tmp_tCC, 'Growth')
        tmp_CC = self.canopy_cover_development(self.var.CC0, self.var.CCx, self.var.CGC, self.var.CDC, tmp_tCC, 'Growth')
        self.var.CC[cond522111] = tmp_CC[cond522111]

        # Determine time required to reach CC on previous day, given CGCadj
        # value
        cond522112 = (cond52211 & np.logical_not(cond522111))
        tReq = self.canopy_cover_required_time(self.var.CC0adj, CCxAdj, CGCadj, self.var.CDC, dtCC, tCCadj, 'CGC')
        tmp_tCC = tReq + dtCC

        # Determine new canopy size
        cond5221121 = (cond522112 & (tmp_tCC > 0))
        # tmp_CC = self.canopy_cover_development(tmp_tCC, 'Growth')
        tmp_CC = self.canopy_cover_development(self.var.CC0adj, CCxAdj, CGCadj, self.var.CDC, tmp_tCC, 'Growth')
        self.var.CC[cond5221121] = tmp_CC[cond5221121]

        # No canopy growth (line 110)
        cond5221122 = (cond522112 & np.logical_not(cond5221121))
        self.var.CC[cond5221122] = self.var.CCprev[cond5221122]

        # No canopy growth (line 115)
        cond52212 = (cond5221 & np.logical_not(cond52211))
        self.var.CC[cond52212] = self.var.CCprev[cond52212]

        # No canopy growth (line 119)
        cond5222 = (cond522 & np.logical_not(cond5221))
        self.var.CC[cond5222] = self.var.CCprev[cond5222]

        # Update CC0 if current canopy cover if less than initial canopy cover size at planting
        cond52221 = (cond5222 & (self.var.CC < self.var.CC0adj))
        self.var.CC0adj[cond52221] = self.var.CC[cond52221]

        # Update actual maximum canopy cover size during growing season
        cond53 = (cond5 & (self.var.CC > self.var.CCxAct))
        self.var.CCxAct[cond53] = self.var.CC[cond53]

        # No more canopy growth is possible or canopy is in decline (line 132)
        cond6 = (self.var.GrowingSeasonIndex & np.logical_not(cond4 | cond5) & (tCCadj > self.var.CanopyDevEnd))

        # Mid-season stage - no canopy growth: update actual maximum canopy
        # cover size during growing season only (i.e. do not update CC)
        cond61 = (cond6 & (tCCadj < self.var.Senescence))
        self.var.CC[cond61] = self.var.CCprev[cond61]
        cond611 = (cond61 & (self.var.CC > self.var.CCxAct))
        self.var.CCxAct[cond611] = self.var.CC[cond611]

        # Late season stage - canopy decline: update canopy decline coefficient
        # for difference between actual and potential CCx, and determine new
        # canopy size
        cond62 = (cond6 & np.logical_not(cond61))
        CDCadj = (self.var.CDC * np.divide(self.var.CCxAct, self.var.CCx, out=np.zeros_like(self.var.CCx), where=self.var.CCx!=0))
        tmp_tCC = tCCadj - self.var.Senescence
        tmp_CC = self.canopy_cover_development(self.var.CC0adj, self.var.CCxAct, self.var.CGC, CDCadj, tmp_tCC, 'Decline')
        self.var.CC[cond62] = tmp_CC[cond62]

    def update_CC_after_senescence(self, tCCadj, dtCC):

        CCsen = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        cond7 = (self.var.GrowingSeasonIndex & (tCCadj >= self.var.Emergence))
        
        # Check for early canopy senescence starting/continuing due to severe
        # water stress
        cond71 = (cond7 & ((tCCadj < self.var.Senescence) | (self.var.tEarlySen > 0)))
        
        # Early canopy senescence
        cond711 = (cond71 & (self.var.Ksw_Sen < 1))

        CDCadj = np.zeros_like(self.var.CDC)
        cond7112 = (cond711 & (self.var.Ksw_Sen > 0.99999))
        CDCadj[cond7112] = 0.0001
        cond7113 = (cond711 & np.logical_not(cond7112))

        CDCadj[cond7113] = ((1. - (self.var.Ksw_Sen ** 8)) * self.var.CDC)[cond7113]

        # Get new canopy cover size after senescence
        cond7114 = (cond711 & (self.var.CCxEarlySen < 0.001))
        CCsen[cond7114] = 0

        # Get time required to reach CC at end of previous day, given CDCadj
        cond7115 = (cond711 & np.logical_not(cond7114))
        tReq = self.canopy_cover_required_time(self.var.CC0adj, self.var.CCxEarlySen, self.var.CGC, CDCadj, dtCC, tCCadj, 'CDC')

        # Calculate GDD's for canopy decline and determine new canopy size
        tmp_tCC = tReq + dtCC
        tmp_CCsen = self.canopy_cover_development(
            self.var.CC0adj,
            self.var.CCxEarlySen,
            self.var.CGC,
            CDCadj, tmp_tCC, 'Decline')
        CCsen[cond7115] = tmp_CCsen[cond7115]

        # Update canopy cover size
        cond7116 = (cond711 & (tCCadj < self.var.Senescence))

        # Limit CC to CCx
        CCsen[cond7116] = np.clip(CCsen, None, self.var.CCx)[cond7116]

        # CC cannot be greater than value on previous day
        self.var.CC[cond7116] = CCsen[cond7116]
        self.var.CC[cond7116] = np.clip(self.var.CC, None, self.var.CCprev)[cond7116]

        # Update maximum canopy cover size during growing season
        self.var.CCxAct[cond7116] = self.var.CC[cond7116]

        # Update CC0 if current CC is less than initial canopy cover size at
        # planting
        cond71161 = (cond7116 & (self.var.CC < self.var.CC0))
        self.var.CC0adj[cond71161] = self.var.CC[cond71161]
        cond71162 = (cond7116 & np.logical_not(cond71161))
        self.var.CC0adj[cond71162] = self.var.CC0[cond71162]

        # Update CC to account for canopy cover senescence due to water stress
        cond7117 = (cond711 & np.logical_not(cond7116))
        self.var.CC[cond7117] = np.clip(self.var.CC, None, CCsen)[cond7117]

        # Check for crop growth termination
        cond7118 = (cond711 & ((self.var.CC < 0.001) & np.logical_not(self.var.CropDead)))
        self.var.CC[cond7118] = 0
        self.var.CropDead[cond7118] = True

        # Otherwise there is no water stress
        cond712 = (cond71 & np.logical_not(cond711))
        self.var.PrematSenes[cond712] = False

        # Rewatering of canopy in late season: get adjusted values of CCx and
        # CDC and update CC
        cond7121 = (cond712 & ((tCCadj > self.var.Senescence) & (self.var.tEarlySen > 0)))
        tmp_tCC = tCCadj - dtCC - self.var.Senescence
        CCxAdj,CDCadj = self.update_CCx_and_CDC(tmp_tCC)
        tmp_tCC = tCCadj - self.var.Senescence
        tmp_CC = self.canopy_cover_development(self.var.CC0adj, CCxAdj, self.var.CGC, CDCadj, tmp_tCC, 'Decline')
        self.var.CC[cond7121] = tmp_CC[cond7121]

        # Check for crop growth termination
        cond71211 = (cond7121 & ((self.var.CC < 0.001) & np.logical_not(self.var.CropDead)))
        self.var.CC[cond71211] = 0
        self.var.CropDead[cond71211] = True

        # Reset early senescence counter
        self.var.tEarlySen[cond712] = 0

        # Adjust CCx for effects of withered canopy
        self.var.CCxW[cond71] = np.clip(self.var.CCxW, self.var.CC, None)[cond71]
        # return CC, CC0adj, CCxAct, CCxW, CropDead, PrematSenes, tEarlySen

    def update_CCx_and_CDC(self, dt):
        """Function to update CCx and CDC parameter values for 
        rewatering in late season of an early declining canopy
        """
        CCxAdj = self.var.CCprev / (1 - 0.05 * (np.exp(dt * (np.divide(self.var.CDC, self.var.CCx, out=np.zeros_like(self.var.CCx), where=self.var.CCx!=0))) - 1))
        CDCadj = self.var.CDC * np.divide(CCxAdj, self.var.CCx, out=np.zeros_like(self.var.CCx), where=self.var.CCx!=0)
        return CCxAdj,CDCadj
        
    def dynamic(self):
        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()            

        # Preallocate some variables
        CCxAdj = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        CDCadj = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        CCsen = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

        # Store initial condition
        self.var.CCprev = np.copy(self.var.CC)

        # Get canopy cover growth over time
        if self.var.CalendarType == 1:
            tCC = np.copy(self.var.DAP)
            dtCC = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
            tCCadj = self.var.DAP - self.var.DelayedCDs
        elif self.var.CalendarType == 2:
            tCC = np.copy(self.var.GDDcum)
            dtCC = np.copy(self.var.GDD)
            tCCadj = self.var.GDDcum - self.var.DelayedGDDs

        # Canopy development (potential)
        self.potential_canopy_development(tCC, dtCC)
        
        # Canopy development (actual)
        self.actual_canopy_development(tCC, tCCadj, dtCC)
        
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

        # **TEST**
        self.var.water_stress_module.dynamic(beta=False)
        
        self.update_CC_after_senescence(tCCadj, dtCC)

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

        self.var.CC[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCadj[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CC_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCadj_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCxW[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCxAct[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCxW_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.CCxAct_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        
