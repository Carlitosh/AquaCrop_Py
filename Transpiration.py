#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class Transpiration(object):
    def __init__(self, Transpiration_variable):
        self.var = Transpiration_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        arr_ones = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.Ksa_Aer = np.copy(arr_zeros)
        self.var.TrPot0 = np.copy(arr_zeros)
        self.var.TrPot_NS = np.copy(arr_zeros)
        self.var.TrAct = np.copy(arr_zeros)
        self.var.TrAct0 = np.copy(arr_zeros)
        self.var.AerDays = np.copy(arr_zeros)
        self.var.AerDaysComp  = np.zeros((self.var.nCrop, self.var.nComp, self.var.nLat, self.var.nLon))
        self.var.Tpot = np.copy(arr_zeros)        
        self.var.TrRatio = np.copy(arr_ones)
        self.var.DaySubmerged = np.copy(arr_zeros)

    def reset_initial_conditions(self):
        cond = self.var.GrowingSeasonDayOne
        cond_comp = np.broadcast_to(cond, self.var.AerDaysComp.shape)
        self.var.AerDays[cond] = 0
        self.var.AerDaysComp[cond_comp] = 0        
        self.var.Tpot[cond] = 0
        self.var.TrRatio[cond] = 1
        # self.var.TrAct[cond] = 0  # TEMP - may not require?
        self.var.DaySubmerged[cond] = 0
        
    def aeration_stress(self):
        """Function to calculate aeration stress coefficient"""

        # Determine aeration stress (root zone)
        cond1 = (self.var.thRZ_Act > self.var.thRZ_Aer)

        # Calculate aeration stress coefficient
        self.var.Ksa_Aer = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
        cond11 = (cond1 & (self.var.AerDays < self.var.LagAer))
        x1 = (self.var.thRZ_Sat - self.var.thRZ_Act)
        x2 = (self.var.thRZ_Sat - self.var.thRZ_Aer)
        stress = (1 - np.divide(x1, x2, out=np.zeros_like(x1), where=x2!=0))
        self.var.Ksa_Aer[cond11] = (1 - ((self.var.AerDays / 3) * stress))[cond11]
        cond12 = (cond1 & np.logical_not(cond11))
        self.var.Ksa_Aer[cond12] = np.divide(x1, x2, out=np.zeros_like(x2), where=x2!=0)[cond12]

        # Increment aeration days counter, or set to zero if there is no stress
        self.var.AerDays[cond1] += 1
        self.var.AerDays[np.logical_not(cond1)] = 0
        self.var.AerDays = np.clip(self.var.AerDays, None, self.var.LagAer)

    def surface_transpiration(self, TrPot):
        
        cond10 = (self.var.GrowingSeasonIndex & (self.var.SurfaceStorage > 0) & (self.var.DaySubmerged < self.var.LagAer))
        self.var.DaySubmerged[cond10] += 1

        # Update anaerobic conditions counter for each compartment
        cond10_comp = np.broadcast_to(cond10, self.var.AerDaysComp.shape)
        self.var.AerDaysComp[cond10_comp] += 1 
        self.var.LagAer_comp = np.broadcast_to(self.var.LagAer, self.var.AerDaysComp.shape)
        self.var.AerDaysComp[cond10_comp] = np.clip(self.var.AerDaysComp, None, self.var.LagAer_comp)[cond10_comp]

        # Reduce actual transpiration that is possible to account for aeration
        # stress due to extended submergence
        fSub = 1 - np.divide(self.var.DaySubmerged.astype(np.float64), self.var.LagAer.astype(np.float64), out=np.zeros_like(self.var.LagAer.astype(np.float64)), where=self.var.LagAer!=0)

        # Transpiration occurs from surface storage
        cond101 = (cond10 & (self.var.SurfaceStorage > (fSub * self.var.TrPot0)))
        self.var.SurfaceStorage[cond101] -= (fSub * self.var.TrPot0)[cond101]
        self.var.TrAct0[cond101] = (fSub * self.var.TrPot0)[cond101]

        # Otherwise there is no transpiration from surface storage
        cond102 = (cond10 & np.logical_not(cond101))
        self.var.TrAct0[cond102] = 0

        # More water can be extracted from soil profile for transpiration
        cond103 = (cond10 & (self.var.TrAct0 < (fSub * self.var.TrPot0)))
        TrPot[cond103] = ((fSub * self.var.TrPot0) - self.var.TrAct0)[cond103]

        # Otherwise no more transpiration possible on current day
        cond104 = (cond10 & np.logical_not(cond103))
        TrPot[cond104] = 0

        cond11 = (self.var.GrowingSeasonIndex & np.logical_not(cond10))
        TrPot[cond11] = self.var.TrPot0[cond11]
        self.var.TrAct0[cond11] = 0
        
    def dynamic(self):
        """Function to calculate crop transpiration on current day"""

        # reset initial conditions
        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()

        # Add crop dimension to ET0
        et0 = self.var.referencePotET[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
        
        # potential transpiration
        # #######################
        
        # 1. No prior water stress
        Kcb = np.copy(self.var.Kcb)
        cond2 = (self.var.GrowingSeasonIndex & (self.var.AgeDays_NS > 5))
        Kcb[cond2] = (self.var.Kcb - ((self.var.AgeDays_NS - 5.) * (self.var.fage / 100.)) * self.var.CCxW_NS)[cond2]
        
        conc = np.copy(self.var.CurrentConc)
        # conc = self.var.conc[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
        co2_mult = (1. - 0.05 * ((conc - self.var.RefConc) / (550. - self.var.RefConc)))
        
        cond4 = (self.var.GrowingSeasonIndex & (conc > self.var.RefConc))
        Kcb[cond4] *= co2_mult[cond4]
        self.var.TrPot_NS = Kcb * self.var.CCadj_NS * et0 * self.var.GrowingSeasonIndex
        
        # 2. Potential prior water stress and/or delayed development
        Kcb = np.copy(self.var.Kcb)
        cond2 = (self.var.GrowingSeasonIndex & (self.var.AgeDays > 5))
        Kcb[cond2] = (self.var.Kcb - ((self.var.AgeDays - 5.) * (self.var.fage / 100)) * self.var.CCxW)[cond2]
        Kcb[cond4] *= co2_mult[cond4]
        self.var.TrPot0 = Kcb * self.var.CCadj * et0 * self.var.GrowingSeasonIndex
        # print 'TrPot_NS   %f' % self.var.TrPot_NS[0,0,0] # OK
        
        # Correct potential transpiration for dying green canopy effects
        cond9 = (self.var.GrowingSeasonIndex & (self.var.CC < self.var.CCxW))
        cond91 = (cond9 & (self.var.CCxW > 0.001) & (self.var.CC > 0.001))
        x = np.divide(self.var.CC, self.var.CCxW, out=np.zeros_like(self.var.CCxW), where=self.var.CCxW!=0)
        self.var.TrPot0[cond91] *= (x ** self.var.a_Tr)[cond91]

        # print 'Kcb        %f' % Kcb[0,0,0]
        # print 'CCadj      %f' % self.var.CCadj[0,0,0]
        # print 'TrPot0     %f' % self.var.TrPot0[0,0,0]

        # surface layer transpiration
        # ###########################
        
        # Potential transpiration counter
        TrPot = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.surface_transpiration(TrPot)
        
        # Update potential root zone transpiration for water stress
        # Maximum stress effect
        self.aeration_stress()
        Ks = np.minimum(self.var.Ksw_StoLin, self.var.Ksa_Aer)
        
        # Update potential transpiration in root zone
        cond11 = (self.var.GrowingSeasonIndex & (self.var.IrrMethod != 4))
        TrPot[cond11] *= Ks[cond11]

        # Maximum sink term
        # #################

        SxComp = np.zeros((self.var.nCrop, self.var.nComp, self.var.nLat, self.var.nLon))

        rootdepth = np.maximum(self.var.Zmin, self.var.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        rootdepth_comp = np.broadcast_to(rootdepth[:,None,:,:], self.var.th.shape)
        # dz = self.var.dz[:,None,None,None] * np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))[None,:,:,:]
        # dzsum = self.var.dzsum[:,None,None,None] * np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))[None,:,:,:]
        comp_sto = (np.round((self.var.dzsum_xy - self.var.dz_xy) * 1000) < np.round(rootdepth_comp * 1000))

        # Fraction of compartment covered by root zone (zero in compartments
        # NOT covered by the root zone)
        RootFact = 1 - ((self.var.dzsum_xy - rootdepth_comp) / self.var.dz_xy)
        RootFact = np.clip(RootFact, 0, 1) * comp_sto

        # Net irrigation mode
        cond12 = (self.var.GrowingSeasonIndex & (self.var.IrrMethod == 4))
        cond12_comp = np.broadcast_to(cond12, SxComp.shape)
        SxComp[cond12_comp] = np.broadcast_to(((self.var.SxTop + self.var.SxBot) / 2.), SxComp.shape)[cond12_comp]

        # Otherwise sink term declines linearly with depth
        cond13 = (self.var.GrowingSeasonIndex & np.logical_not(cond12))

        if (np.any(cond13)):
            comp = 0
            comp_sto_sum = np.sum(comp_sto, axis=1)
            SxCompBot = np.copy(self.var.SxTop)
            while np.any(comp < comp_sto_sum):
                SxCompTop = np.copy(SxCompBot)
                cond131 = (cond13 & (self.var.dzsum_xy[:,comp,...] <= rootdepth))
                SxCompBot[cond131] = (self.var.SxBot * self.var.rCor + ((self.var.SxTop - self.var.SxBot * self.var.rCor) * ((rootdepth - self.var.dzsum_xy[:,comp,...]) / rootdepth)))[cond131]
                cond132 = (cond13 & np.logical_not(cond131))
                SxCompBot[cond132] = (self.var.SxBot * self.var.rCor)[cond132]
                SxComp[:,comp,...][cond13] = ((SxCompTop + SxCompBot) / 2)[cond13]
                comp += 1

        SxComp *= comp_sto

        # Extract water
        self.var.TrAct.fill(0.)
        # print 'TrPot   %f' % TrPot[0,0,0]
        ToExtract = np.copy(TrPot)
        cond14_ini = (self.var.GrowingSeasonIndex & (ToExtract > 0))
        if (np.any(cond14_ini)):
            comp = 0
            comp_sto_sum = np.sum(comp_sto, axis=1)
            while np.any((comp < comp_sto_sum) & (ToExtract > 0)):

                cond14 = (self.var.GrowingSeasonIndex & (comp_sto[:,comp,...]) & (ToExtract > 0))

                # Determine TAW for compartment
                thTAW = self.var.th_fc_comp[:,comp,...] - self.var.th_wp_comp[:,comp,...]
                p_up_sto = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
                cond141 = (cond14 & (self.var.ETadj == 1))
                p_up_sto[cond141] = (self.var.p_up2 + (0.04 * (5. - et0)) * (np.log10(10. - 9. * self.var.p_up2)))[cond141]

                # Determine critical water content at which stomatal closure
                # will occur in compartment
                thCrit = (self.var.th_fc_comp[:,comp,...] - (thTAW * p_up_sto))

                # Check for soil water stress
                KsComp = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

                # No stress
                cond142 = (cond14 & (self.var.th[:,comp,...] >= thCrit))
                KsComp[cond142] = 1.

                # Transpiration from compartment is affected by water stress
                cond143 = (cond14 & (self.var.th[:,comp,...] > self.var.th_wp_comp[:,comp,...]) & np.logical_not(cond142))
                Wrel = ((self.var.th_fc_comp[:,comp,...] - self.var.th[:,comp,...]) / (self.var.th_fc_comp[:,comp,...] - self.var.th_wp_comp[:,comp,...]))
                pRel = ((Wrel - self.var.p_up2) / (self.var.p_lo2 - self.var.p_up2))
                KsComp[cond143] = (1 - ((np.exp(pRel * self.var.fshape_w2) - 1) / (np.exp(self.var.fshape_w2) - 1)))[cond143]
                KsComp = np.clip(KsComp, 0, 1)
                KsComp[pRel <= 0] = 1
                KsComp[pRel >= 1] = 0

                # Otherwise no transpiration is possible from the compartment
                # as water does not exceed wilting point
                KsComp[(cond14 & np.logical_not(cond142 | cond143))] = 0

                # Adjust compartment stress factor for aeration stress
                AerComp = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

                # Full aeration stress - no transpiration possible from
                # compartment
                cond144 = (cond14 & (self.var.DaySubmerged >= self.var.LagAer))
                cond145 = (cond14 & (self.var.th[:,comp,...] > (self.var.th_s_comp[:,comp,...] - (self.var.Aer / 100))) & np.logical_not(cond144))
                self.var.AerDaysComp[:,comp,...][cond145] += 1
                fAer = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
                cond1451 = (cond145 & (self.var.AerDaysComp[:,comp,...] >= self.var.LagAer))
                self.var.AerDaysComp[:,comp,...][cond1451] = self.var.LagAer[cond1451]
                fAer[cond1451] = 0

                # Calculate aeration stress factor
                AerComp[cond145] = ((self.var.th_s_comp[:,comp,...] - self.var.th[:,comp,...]) / (self.var.th_s_comp[:,comp,...] - (self.var.th_s_comp[:,comp,...] - (self.var.Aer / 100))))[cond145]
                AerComp = np.clip(AerComp, 0, None)
                AerComp_divd = (fAer + (self.var.AerDaysComp[:,comp,...] - 1) * AerComp)
                AerComp_divs = (fAer + self.var.AerDaysComp[:,comp,...] - 1)
                AerComp[cond145] = np.divide(AerComp_divd, AerComp_divs, out=np.zeros_like(AerComp_divs), where=AerComp_divs!=0)[cond145]

                # Otherwise there is no aeration stress as number of submerged
                # days does not exceed threshold for initiation of aeration
                # stress
                cond146 = (cond14 & np.logical_not(cond144 | cond145))
                AerComp[cond146] = 1
                self.var.AerDaysComp[:,comp,...][cond146] = 0

                # Extract water
                ThToExtract = ((ToExtract / 1000) / self.var.dz[comp])
                Sink = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

                # Don't reduce compartment sink for stomatal water stress if in
                # net irrigation mode. Stress only occurs due to deficient
                # aeration conditions
                cond147 = (cond14 & self.var.IrrMethod == 4)
                Sink[cond147] = (AerComp * SxComp[:,comp,...] * RootFact[:,comp,...])[cond147]

                # Otherwise, reduce compartment sink for greatest of stomatal
                # and aeration stress
                cond148 = (cond14 & np.logical_not(cond147))
                cond1481 = (cond148 & (KsComp == AerComp))
                Sink[cond1481] = (KsComp * SxComp[:,comp,...] * RootFact[:,comp,...])[cond1481]
                cond1482 = (cond148 & np.logical_not(cond1481))
                Sink[cond1482] = (np.minimum(KsComp,AerComp) * SxComp[:,comp,...] * RootFact[:,comp,...])[cond1482]

                # Limit extraction to demand
                Sink = np.clip(Sink, None, ThToExtract)

                # Limit extraction to avoid compartment water content dropping
                # below air dry
                cond149 = (cond14 & ((self.var.th[:,comp,...] - Sink) < self.var.th_dry_comp[:,comp,...]))
                Sink[cond149] = (self.var.th[:,comp,...] - self.var.th_dry_comp[:,comp,...])[cond149]
                Sink = np.clip(Sink, 0, None)

                # Update water content in compartment
                self.var.th[:,comp,...][cond14] -= Sink[cond14]

                # Update amount of water to extract
                ToExtract[cond14] -= (Sink * 1000 * self.var.dz[comp])[cond14]

                # Update actual transpiration
                # TrActComp[comp,:][cond14] += (Sink * 1000 * dz[comp])[cond14]
                self.var.TrAct[cond14] += (Sink * 1000 * self.var.dz[comp])[cond14]

                # Update compartment counter
                comp += 1

        # Add net irrigation water requirement
        # ####################################
        self.var.root_zone_water_module.dynamic()
        cond15 = (self.var.GrowingSeasonIndex & (self.var.IrrMethod == 4) & (TrPot > 0))
        thCrit = self.var.thRZ_Wp + ((self.var.NetIrrSMT / 100) * (self.var.thRZ_Fc - self.var.thRZ_Wp))
        cond151 = (cond15 & (self.var.thRZ_Act < thCrit))
        cond151_comp = np.broadcast_to(cond151, self.var.th.shape)

        # Calculate thCrit in each compartment
        thCrit_comp = (self.var.th_wp_comp + ((self.var.NetIrrSMT / 100) * (self.var.th_fc_comp - self.var.th_wp_comp)))

        # Determine necessary change in water content in compartments to reach
        # critical water content
        dWC = RootFact * (thCrit_comp - self.var.th * 1000 * self.var.dz_xy)
        self.var.th[cond151_comp] = (self.var.th + (dWC / (1000 * self.var.dz_xy)))[cond151_comp]
        self.var.IrrNet[cond151] = np.sum(dWC, axis=1)[cond151]
        
        # Update net irrigation counter for the growing season
        self.var.IrrNetCum += self.var.IrrNet
        
        # Add any surface transpiration to root zone total
        self.var.TrAct += self.var.TrAct0

        # Update transpiration ratio
        # ##########################        
        cond17 = (self.var.GrowingSeasonIndex & (self.var.TrPot0 > 0))
        cond171 = (cond17 & (self.var.TrAct < self.var.TrPot0))
        self.var.TrRatio[cond171] = np.divide(self.var.TrAct, self.var.TrPot0, out=np.zeros_like(self.var.TrPot0), where=self.var.TrPot0!=0)[cond171]
        cond172 = (cond17 & np.logical_not(cond171))
        self.var.TrRatio[cond172] = 1
        cond18 = (self.var.GrowingSeasonIndex & np.logical_not(cond17))
        self.var.TrRatio[cond18] = 1
        self.var.TrRatio = np.clip(self.var.TrRatio, 0, 1)

        # No transpiration or irrigation if outside growing season
        # ########################################################
        self.var.TrAct[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.TrPot0[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.TrPot_NS[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.IrrNet[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.IrrNetCum[np.logical_not(self.var.GrowingSeasonIndex)] = 0

        # Store potential transpiration for irrigation calculations on next day
        self.var.Tpot = np.copy(self.var.TrPot0)

        # print 'TrAct   %f' % self.var.TrAct[0,0,0]
