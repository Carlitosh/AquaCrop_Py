#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class SoilEvaporation(object):
    """Class to represent daily soil evaporation"""

    def __init__(self, SoilEvaporation_variable):
        self.var = SoilEvaporation_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.Epot = np.copy(arr_zeros)
        self.var.Stage2 = np.copy(arr_zeros.astype(bool))
        self.var.EvapZ = np.copy(arr_zeros)
        self.var.Wstage2 = np.copy(arr_zeros)
        self.var.Wsurf = np.copy(arr_zeros)
        self.var.Wevap_Act = np.copy(arr_zeros)
        self.var.Wevap_Sat = np.copy(arr_zeros)
        self.var.Wevap_Fc = np.copy(arr_zeros)
        self.var.Wevap_Wp = np.copy(arr_zeros)
        self.var.Wevap_Dry = np.copy(arr_zeros)

    def reset_initial_conditions(self):
        pass

    def evap_layer_water_content(self):
        """Function to get water contents in the evaporation layer"""
        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()
        
        arr_ones = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
        dz = self.var.dz[None,:,None,None] * arr_ones[:,None,:,:]
        dzsum = self.var.dzsum[None,:,None,None] * arr_ones[:,None,:,:]

        # Find compartments covered by evaporation layer
        evapz_comp = self.var.EvapZ[:,None,:,:] * np.ones((self.var.nComp))[None,:,None,None]
        comp_sto = (np.round((dzsum - dz) * 1000) < np.round(evapz_comp * 1000))
        factor = 1 - ((dzsum - evapz_comp) / dz)
        factor = np.clip(factor, 0, 1) * comp_sto

        # Water storages in evaporation layer (mm)
        Wevap_Act = np.sum((factor * 1000 * self.var.th * dz), axis=1)
        self.var.Wevap_Act = np.clip(Wevap_Act, 0, None)
        self.var.Wevap_Sat = np.sum((factor * 1000 * self.var.th_s_comp * dz), axis=1)
        self.var.Wevap_Fc = np.sum((factor * 1000 * self.var.th_fc_comp * dz), axis=1)
        self.var.Wevap_Wp = np.sum((factor * 1000 * self.var.th_wp_comp * dz), axis=1)
        self.var.Wevap_Dry = np.sum((factor * 1000 * self.var.th_dry_comp * dz), axis=1)

    def prepare_stage_two_evaporation(self):
        self.evap_layer_water_content()
        Wstage2 = ((self.var.Wevap_Act - (self.var.Wevap_Fc - self.var.REW)) / (self.var.Wevap_Sat - (self.var.Wevap_Fc - self.var.REW)))
        Wstage2 = (np.round((Wstage2 * 100)) / 100)
        Wstage2 = np.clip(Wstage2, 0, None)
        # self.var.Wstage2 = Wstage2
        return Wstage2

    def potential_soil_evaporation_rate(self, tAdj):

        # No canopy cover outside of growing season so potential soil
        # evaporation only depends on reference evapotranspiration
        et0 = self.var.referencePotET[None,:,:] * np.ones((self.var.nCrop))[:,None,None]        
        EsPot = (self.var.Kex * et0)

        # Calculate maximum potential soil evaporation and potential soil
        # evaporation given current canopy size
        EsPotMax = (self.var.Kex * et0 * (1 - self.var.CCxW * (self.var.fwcc / 100)))
        EsPot[self.var.GrowingSeasonIndex] = (self.var.Kex * (1 - self.var.CCadj) * et0)[self.var.GrowingSeasonIndex]

        # Adjust potential soil evaporation for effects of withered canopy
        cond3 = (self.var.GrowingSeasonIndex & (tAdj > self.var.Senescence) & (self.var.CCxAct > 0))
        mult = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
        cond31 = (cond3 & (self.var.CC > (self.var.CCxAct / 2)))
        cond311 = (cond31 & (self.var.CC > self.var.CCxAct))
        mult[cond311] = 0
        cond312 = (cond31 & np.logical_not(cond311))
        mult_divd = (self.var.CCxAct - self.var.CC)
        mult_divs = (self.var.CCxAct / 2)
        mult[cond312] = np.divide(mult_divd, mult_divs, out=np.zeros_like(mult_divs), where=mult_divs!=0)[cond312]

        EsPot[cond3] = (EsPot * (1 - self.var.CCxAct * (self.var.fwcc / 100) * mult))[cond3]
        CCxActAdj = ((1.72 * self.var.CCxAct) + (self.var.CCxAct ** 2) - 0.3 * (self.var.CCxAct ** 3))
        EsPotMin = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        EsPotMin[cond3] = (self.var.Kex * (1 - CCxActAdj) * et0)[cond3]
        EsPotMin = np.clip(EsPotMin, 0, None)

        # Line 85-89 of AOS_SoilEvaporation.m
        EsPot[cond3] = np.clip(EsPot, EsPotMin, EsPotMax)[cond3]

        cond4 = (self.var.GrowingSeasonIndex & self.var.PrematSenes)
        EsPot[cond4] = np.clip(EsPot, None, EsPotMax)[cond4]
        # self.EsPot = EsPot
        return EsPot

    def adjust_potential_soil_evaporation_for_irrigation(self, EsPot):
        EsPotIrr = np.copy(EsPot)
        cond1 = ((irr_depth > 0) & (irr_method != 4))
        cond11 = (cond1 & ((prec > 0) | (surface_storage > 0)))
        EsPotIrr[cond11] = EsPot[cond11]
        cond12 = (cond1 & np.logical_not(cond11))
        EsPotIrr[cond12] = (EsPot * (wet_surf / 100))[cond12]  # TODO: more informative name for wet_surf
        return EsPotIrr

    def adjust_potential_soil_evaporation_for_mulches(self, EsPot):

        # NB if surface is flooded then there is no adjustment of potential soil
        # evaporation, regardless of mulches    
        EsPotMul = np.copy(EsPot)
        cond1 = (self.var.SurfaceStorage < 0.000001)
        cond11 = (cond1 & (self.var.Mulches == 1))
        cond111 = (cond11 & self.var.GrowingSeasonIndex)
        EsPotMul[cond111] = (EsPot * (1 - self.var.fMulch * (self.var.MulchPctGS / 100)))[cond111]
        cond112 = (cond11 & np.logical_not(self.var.GrowingSeasonIndex))
        EsPotMul[cond112] = (EsPot * (1 - self.var.fMulch * (self.var.MulchPctOS / 100)))[cond112]
        return EsPotMul

    def extract_water(self, ToExtract, ToExtractStg):
        arr_ones = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        dz = self.var.dz[None,:,None,None] * arr_ones[:,None,:,:]
        dzsum = self.var.dzsum[None,:,None,None] * arr_ones[:,None,:,:]

        # Determine fraction of compartments covered by evaporation layer
        evapzmin_comp = self.var.EvapZmin[:,None,:,:] * np.ones((self.var.nComp))[None,:,None,None]
        comp_sto = (np.round((dzsum - dz) * 1000) < np.round(evapzmin_comp * 1000))
        factor = 1 - ((dzsum - evapzmin_comp) / dz)
        factor = np.clip(factor, 0, 1) * comp_sto

        comp_sto = np.sum(comp_sto, axis=1)
        comp = 0
        while np.any((comp < comp_sto) & (ToExtractStg > 0) & (ToExtract > 0)):

            cond101 = ((comp < comp_sto) & (ToExtractStg > 0) & (ToExtract > 0))

            # Water available in compartment for extraction (mm)
            Wdry = 1000 * self.var.th_dry_comp[:,comp,...] * dz[:,comp,...]  
            W = 1000 * self.var.th[:,comp,...] * dz[:,comp,...]
            AvW = np.copy(arr_zeros)
            AvW[cond101] = ((W - Wdry) * factor[:,comp,...])[cond101]
            AvW = np.clip(AvW, 0, None)

            # Determine amount by which to adjust variables
            cond1011 = (cond101 & (AvW >= ToExtractStg))
            self.var.EsAct[cond1011] += ToExtractStg[cond1011]
            W[cond1011] -= ToExtractStg[cond1011]
            ToExtract[cond1011] -= ToExtractStg[cond1011]
            ToExtractStg[cond1011] = 0

            cond1012 = (cond101 & np.logical_not(cond1011))
            self.var.EsAct[cond1012] += AvW[cond1012]
            W[cond1012] -= AvW[cond1012]
            ToExtract[cond1012] -= AvW[cond1012]
            ToExtractStg[cond1012] -= AvW[cond1012]

            # Update water content
            self.var.th[:,comp,...][cond101] = (W / (1000 * dz[:,comp,...]))[cond101]
            comp += 1

    def relative_depletion(self):

        # Get current water storage
        self.evap_layer_water_content()

        # Get water storage (mm) at start of stage 2 evaporation
        Wupper = (self.var.Wstage2 * (self.var.Wevap_Sat - (self.var.Wevap_Fc - self.var.REW)) + (self.var.Wevap_Fc - self.var.REW))
        # Get water storage (mm) when there is no evaporation
        Wlower = np.copy(self.var.Wevap_Dry)
        # Get relative depletion of evaporation storage in stage 2
        Wrel_divd = (self.var.Wevap_Act - Wlower)
        Wrel_divs = (Wupper - Wlower)
        Wrel = np.divide(Wrel_divd, Wrel_divs, out=np.zeros_like(Wrel_divs), where=Wrel_divs!=0)
        return Wrel, Wlower, Wupper    
            
    def dynamic(self):
        
        # Add crop dimension to self.var.vars
        et0 = self.var.referencePotET[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
        prec = self.var.precipitation[None,:,:] * np.ones((self.var.nCrop))[:,None,None]

        # Prepare stage 2 evaporation (REW gone), if day one of simulation
        cond1 = (self.var._modelTime.timeStepPCR == 1)
        if self.var._modelTime.timeStepPCR == 1:
            self.var.Wsurf.fill(0)
            self.var.EvapZ = np.copy(self.var.EvapZmin)
            self.var.Wstage2 = self.prepare_stage_two_evaporation()
            
        # Prepare soil evaporation stage 1 - adjust water in surface evaporation
        # layer for any infiltration (only do this if rainfall occurs or when
        # irrigation is triggered)
        cond2 = ((prec > 0) | np.any(((self.var.Irr > 0) & (self.var.IrrMethod != 4)), axis=0))
        cond21 = (cond2 & (self.var.Infl > 0))
        self.var.Wsurf[cond21] = self.var.Infl[cond21]
        self.var.Wsurf[cond21] = np.clip(self.var.Wsurf, None, self.var.REW)[cond21]
        self.var.Wstage2[cond21] = 0  # TODO: is this right?
        self.var.EvapZ[cond21] = self.var.EvapZmin[cond21]
        self.var.Stage2[cond21] = False

        # Calculate potential soil evaporation rate (mm/day)
        if self.var.CalendarType == 1:
            tAdj = (self.var.DAP - self.var.DelayedCDs) # * growing_season_index
        elif self.var.CalendarType == 2:
            tAdj = (self.var.GDDcum - self.var.DelayedGDDs) # * growing_season_index

        EsPot = self.potential_soil_evaporation_rate(tAdj)
        
        # Adjust potential soil evaporation for mulches and/or partial wetting
        EsPotMul = np.copy(EsPot)
        cond1 = (self.var.SurfaceStorage < 0.000001)
        cond11 = (cond1 & (self.var.Mulches == 1))
        cond111 = (cond11 & self.var.GrowingSeasonIndex)
        EsPotMul[cond111] = (EsPot * (1 - self.var.fMulch * (self.var.MulchPctGS / 100)))[cond111]
        cond112 = (cond11 & np.logical_not(self.var.GrowingSeasonIndex))
        EsPotMul[cond112] = (EsPot * (1 - self.var.fMulch * (self.var.MulchPctOS / 100)))[cond112]

        # Partial surface wetting by irrigation
        EsPotIrr = np.copy(EsPot)
        cond1 = ((self.var.Irr > 0) & (self.var.IrrMethod != 4))
        cond11 = (cond1 & ((prec > 0) | (self.var.SurfaceStorage > 0)))
        EsPotIrr[cond11] = EsPot[cond11]
        cond12 = (cond1 & np.logical_not(cond11))
        EsPotIrr[cond12] = (EsPot * (self.var.WetSurf / 100))[cond12]  # TODO: more informative name for wet_surf
        # return EsPotIrr

        # Assign minimum value (mulches and partial wetting don't combine)
        EsPot = np.minimum(EsPotIrr, EsPotMul)

        # Initialise actual evaporation counter
        self.var.EsAct = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

        # Surface evaporation        
        # EsActSurf = surface_evaporation(self.var.SurfaceStorage, EsPot)
        EsActSurf = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        cond9 = (self.var.SurfaceStorage > 0)
        cond91 = (cond9 & (self.var.SurfaceStorage > EsPot))
        EsActSurf[cond91] = EsPot[cond91]
        self.var.SurfaceStorage[cond91] -= EsActSurf[cond91]
        cond92 = (cond9 & np.logical_not(cond91))
        EsActSurf[cond92] = self.var.SurfaceStorage[cond92]        
        
        self.var.EsAct += EsActSurf
        cond = ((self.var.SurfaceStorage > 0) & (self.var.SurfaceStorage <= EsPot))
        self.var.SurfaceStorage -= EsActSurf
        self.var.Wsurf[cond] = self.var.REW[cond]
        self.var.Wstage2[cond] = 0
        self.var.EvapZ[cond] = self.var.EvapZmin[cond]
        
        # Stage 1 evaporation
        dz = self.var.dz[None,:,None,None] * np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))[:,None,:,:]
        dzsum = self.var.dzsum[None,:,None,None] * np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))[:,None,:,:]

        # Determine total water to be extracted
        # print EsPot[0,0,0]
        # print self.var.EsAct[0,0,0]
        ToExtract = EsPot - self.var.EsAct
        # print ToExtract[0,0,0]

        # Determine total water to be extracted in stage one (limited by
        # surface layer water storage)
        ExtractPotStg1 = np.minimum(ToExtract, self.var.Wsurf)

        # Extract water
        cond10 = (ExtractPotStg1 > 0)
        self.extract_water(ToExtract, ExtractPotStg1)
        
        # Update surface evaporation layer water balance
        self.var.Wsurf[cond10] -= self.var.EsAct[cond10]
        cond102 = (cond10 & ((self.var.Wsurf < 0) | (ExtractPotStg1 > 0.0001)))
        self.var.Wsurf[cond102] = 0

        # If surface storage completely depleted, prepare stage 2 evaporation
        cond103 = (cond10 & (self.var.Wsurf < 0.0001))
        Wstage2_tmp = self.prepare_stage_two_evaporation()
        self.var.Wstage2[cond103] = Wstage2_tmp[cond103]  # TODO: add MASK argument to function?
        
        # Stage 2 evaporation
        # self.var.th, self.var.EsAct = soil_evaporation_stage_two(
        #     self.var.th,
        #     self.var.th_s_comp,
        #     self.var.th_fc_comp,
        #     self.var.th_wp_comp,
        #     self.var.th_dry_comp,
        #     EsPot,
        #     self.var.EsAct,
        #     self.var.dz,
        #     self.var.dzsum,
        #     self.var.EvapZmin,
        #     self.var.EvapZmax,
        #     self.var.REW,
        #     self.var.EvapZ,
        #     self.var.fWrelExp,
        #     self.var.fevap,
        #     self.var.Wstage2,
        #     ToExtract,
        #     EvapTimeSteps=20)

        # print ToExtract[0,0,0]
        # print self.var.EvapZmax[0,0,0]
        # print self.var.EvapZmin[0,0,0]
        # print self.var.EvapZ[0,0,0]
        # print self.var.fevap[0,0,0]
        # print self.var.fWrelExp[0,0,0]

        # Extract water
        EvapTimeSteps = 20
        cond11 = (ToExtract > 0)        
        if np.any(cond11):
            Edt = ToExtract / EvapTimeSteps
            # Loop sub-daily time steps
            for jj in range(EvapTimeSteps):

                # Get water storage (mm) at start of stage 2 evaporation
                Wrel, Wlower, Wupper = self.relative_depletion()

                # Check if need to expand evaporative layer
                cond111 = (cond11 & (self.var.EvapZmax > self.var.EvapZmin))
                Wcheck = (self.var.fWrelExp * ((self.var.EvapZmax - self.var.EvapZ) / (self.var.EvapZmax - self.var.EvapZmin)))
                while np.any(cond111 & (Wrel < Wcheck) & (self.var.EvapZ < self.var.EvapZmax)):
                    cond1111 = (cond111 & (Wrel < Wcheck) & (self.var.EvapZ < self.var.EvapZmax))

                    # Expand evaporation layer by 1mm
                    self.var.EvapZ[cond1111] += 0.001

                    # Recalculate current water storage for new EvapZ
                    Wrel, Wlower, Wupper = self.relative_depletion()
                    Wcheck = (self.var.fWrelExp * ((self.var.EvapZmax - self.var.EvapZ) / (self.var.EvapZmax - self.var.EvapZmin)))

                # Get stage 2 evaporation reduction coefficient
                Kr = ((np.exp(self.var.fevap * Wrel) - 1) / (np.exp(self.var.fevap) - 1))
                Kr = np.clip(Kr, None, 1)

                # Get water to extract (NB Edt is zero in cells which do not
                # need stage 2, so no need for index)
                # print '********'
                # print self.var.Wevap_Act[0,0,0]
                # print self.var.Wevap_Sat[0,0,0]
                # print self.var.Wevap_Fc[0,0,0]
                # print self.var.Wevap_Dry[0,0,0]
                # print self.var.Wstage2[0,0,0]
                # print Wrel[0,0,0]
                # print Kr[0,0,0]
                # print Edt[0,0,0]
                ToExtractStg2 = (Kr * Edt)
                # print ToExtractStg2[0,0,0]
                self.extract_water(ToExtract, ToExtractStg2)
                # print self.var.th[0,0,0,0]
                # thnew, EsAct, ToExtract, ToExtractStg2 = extract_water(
                #     thnew, th_dry, dz, dzsum, EvapZmin, EsAct, ToExtract, ToExtractStg2)

        # print self.var.th[0,0,0,0]
        # return thnew, EsAct

        # Store potential evaporation for irrigation calculations on next day
        self.var.Epot = np.copy(EsPot)

        # print self.var.th[0,0,0,0]
