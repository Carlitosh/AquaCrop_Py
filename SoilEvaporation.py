#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from hydrology_funs import *

import logging
logger = logging.getLogger(__name__)

class SoilEvaporation(object):
    """Class to represent daily soil evaporation"""

    def __init__(self, SoilEvaporation_variable):
        self.var = SoilEvaporation_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.Wevap_Act = np.copy(arr_zeros)  # TODO: create Wevap class
        self.Wevap_Sat = np.copy(arr_zeros)
        self.Wevap_Fc = np.copy(arr_zeros)
        self.Wevap_Wp = np.copy(arr_zeros)
        self.Wevap_Dry = np.copy(arr_zeros)
        
    def dynamic(self):
        
        # Add rotation dimension to self.var.vars
        et0 = self.var.referencePotET[None,:,:] * np.ones((self.var.nRotation))[:,None,None]
        prec = self.var.precipitation[None,:,:] * np.ones((self.var.nRotation))[:,None,None]

        # Prepare stage 2 evaporation (REW gone), if day one of simulation
        cond1 = (self.var._modelTime.timeStepPCR == 1)
        if self.var._modelTime.timeStepPCR == 1:
            self.var.Wsurf.fill(0)
            self.var.EvapZ = np.copy(self.var.EvapZmin)
            self.var.Wstage2 = prepare_stage_two_evaporation(
                self.var.th,
                self.var.th_s_comp,
                self.var.th_fc_comp,
                self.var.th_wp_comp,
                self.var.th_dry_comp,
                self.var.dz,
                self.var.REW,
                self.var.EvapZ)
            
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
        tAdj = adjust_time_for_delayed_development(
            self.var.CalendarType,
            self.var.DAP,
            self.var.DelayedCDs,
            self.var.DelayedGDDs,
            self.var.GDDcum)

        EsPot = potential_soil_evaporation_rate(
            self.var.GrowingSeasonIndex,
            tAdj,
            self.var.CC,
            self.var.CCadj,
            self.var.CCxAct,
            self.var.CCxW,
            self.var.Senescence,
            self.var.PrematSenes,
            et0,
            self.var.Kex,
            self.var.fwcc)
        
        # Adjust potential soil evaporation for mulches and/or partial wetting
        EsPotMul = adjust_potential_soil_evaporation_for_mulches(
            self.var.GrowingSeasonIndex,
            self.var.SurfaceStorage,
            EsPot,
            self.var.Mulches,
            self.var.fMulch,
            self.var.MulchPctGS,
            self.var.MulchPctOS)

        # Partial surface wetting by irrigation
        EsPotIrr = adjust_potential_soil_evaporation_for_irrigation(
            self.var.SurfaceStorage,
            EsPot,
            self.var.IrrMethod,
            self.var.Irr,
            prec,
            self.var.WetSurf)

        # Assign minimum value (mulches and partial wetting don't combine)
        EsPot = np.minimum(EsPotIrr, EsPotMul)

        # Initialise actual evaporation counter
        self.var.EsAct = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))

        # Surface evaporation        
        EsActSurf = surface_evaporation(self.var.SurfaceStorage, EsPot)
        self.var.EsAct += EsActSurf
        cond = ((self.var.SurfaceStorage > 0) & (self.var.SurfaceStorage <= EsPot))
        self.var.SurfaceStorage -= EsActSurf
        self.var.Wsurf[cond] = self.var.REW[cond]
        self.var.Wstage2[cond] = 0
        self.var.EvapZ[cond] = self.var.EvapZmin[cond]
        
        # Stage 1 evaporation
        self.var.th, self.var.Wsurf, ToExtract, self.var.Wstage2 = soil_evaporation_stage_one(
            self.var.th,
            self.var.th_s_comp,
            self.var.th_fc_comp,
            self.var.th_wp_comp,
            self.var.th_dry_comp,
            EsPot,
            self.var.EsAct,
            self.var.dz,
            self.var.dzsum,
            self.var.EvapZmin,
            self.var.REW,
            self.var.EvapZ,
            self.var.Wsurf,
            self.var.Wstage2)

        # Stage 2 evaporation
        self.var.th, self.var.EsAct = soil_evaporation_stage_two(
            self.var.th,
            self.var.th_s_comp,
            self.var.th_fc_comp,
            self.var.th_wp_comp,
            self.var.th_dry_comp,
            EsPot,
            self.var.EsAct,
            self.var.dz,
            self.var.dzsum,
            self.var.EvapZmin,
            self.var.EvapZmax,
            self.var.REW,
            self.var.EvapZ,
            self.var.fWrelExp,
            self.var.fevap,
            self.var.Wstage2,
            ToExtract,
            EvapTimeSteps=20)
                
        # Store potential evaporation for irrigation calculations on next day
        self.var.Epot = np.copy(EsPot)

        # print self.var.th[0,0,0,0]
