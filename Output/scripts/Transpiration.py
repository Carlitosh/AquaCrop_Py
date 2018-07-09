#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from hydrology_funs import *

import logging
logger = logging.getLogger(__name__)

class Transpiration(object):
    def __init__(self, Transpiration_variable):
        self.var = Transpiration_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.Ksa_Aer = np.copy(arr_zeros)
        self.TrPot0 = np.copy(arr_zeros)
        self.TrPot_NS = np.copy(arr_zeros)
        self.TrAct = np.copy(arr_zeros)
        self.TrAct0 = np.copy(arr_zeros)
        
    def aeration_stress(self):
        """Function to calculate aeration stress coefficient"""

        # Determine aeration stress (root zone)
        cond1 = (self.var.thRZ_Act > self.var.thRZ_Aer)

        # Calculate aeration stress coefficient
        self.var.Ksa_Aer = np.ones((self.var.nRotation, self.var.nLat, self.var.nLon))
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

    def dynamic(self):
        """Function to calculate crop transpiration on current day"""

        # Add rotation dimension to ET0
        et0 = self.var.referencePotET[None,:,:] * np.ones((self.var.nRotation))[:,None,None]
        # arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        
        # potential transpiration
        # #######################
        
        # 1. No prior water stress
        self.var.TrPot_NS = potential_transpiration(
            self.var.GrowingSeasonIndex,
            self.var.Kcb,
            self.var.AgeDays_NS, self.var.fage, self.var.CCadj_NS, self.var.CCxW_NS, self.var.conc, self.var.RefConc, et0)

        # 2. Potential prior water stress and/or delayed development
        self.var.TrPot0 = potential_transpiration(
            self.var.GrowingSeasonIndex,
            self.var.Kcb,
            self.var.AgeDays, self.var.fage, self.var.CCadj, self.var.CCxW, self.var.conc, self.var.RefConc, et0)
        
        # Correct potential transpiration for dying green canopy effects
        cond9 = (self.var.GrowingSeasonIndex & (self.var.CC < self.var.CCxW))
        cond91 = (cond9 & (self.var.CCxW > 0.001) & (self.var.CC > 0.001))
        x = np.divide(self.var.CC, self.var.CCxW, out=np.zeros_like(self.var.CCxW), where=self.var.CCxW!=0)
        self.var.TrPot0[cond91] *= (x ** self.var.a_Tr)[cond91]

        # surface layer transpiration
        # ###########################
        
        # Potential transpiration counter
        TrPot = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        TrPot_surf, self.var.TrAct0, self.var.AerDaysComp = surface_transpiration(
            self.var.GrowingSeasonIndex,
            self.var.SurfaceStorage,
            self.var.DaySubmerged,
            self.var.LagAer,
            self.var.AerDaysComp,
            self.var.TrPot0)
        TrPot += TrPot_surf

        # Update potential root zone transpiration for water stress
        # Maximum stress effect
        self.aeration_stress()
        Ks = np.minimum(self.var.Ksw_StoLin, self.var.Ksa_Aer)
        
        # Update potential transpiration in root zone
        cond11 = (self.var.GrowingSeasonIndex & (self.var.IrrMethod != 4))
        TrPot[cond11] *= Ks[cond11]

        # Maximum sink term
        comp_sto, RootFact, SxComp = maximum_sink_term(
            self.var.GrowingSeasonIndex,
            self.var.IrrMethod,
            self.var.SxTop,
            self.var.SxBot,
            self.var.rCor,
            self.var.Zmin, self.var.Zroot, self.var.dz, self.var.dzsum)

        # Extract water
        self.var.th, self.var.AerDaysComp, self.var.TrAct = extract_transpiration_water(
            TrPot, self.var.GrowingSeasonIndex, comp_sto, self.var.DaySubmerged,
            self.var.IrrMethod, self.var.th, self.var.th_s_comp,
            self.var.th_fc_comp, self.var.th_wp_comp,
            self.var.th_dry_comp, self.var.ETadj, self.var.p_up2,
            self.var.p_lo2, self.var.fshape_w2, self.var.LagAer,
            self.var.Aer, self.var.AerDaysComp, et0, self.var.dz,
            SxComp, RootFact)

        # Add net irrigation water requirement
        # ####################################
        self.var.root_zone_water_module.dynamic()  # TODO - how to handle this???
        self.var.th, self.var.IrrNet = add_net_irrigation(
            self.var.GrowingSeasonIndex,
            self.var.IrrMethod,
            self.var.th,
            self.var.th_wp_comp,
            self.var.th_fc_comp,
            self.var.thRZ_Wp, self.var.thRZ_Fc, self.var.thRZ_Act, self.var.NetIrrSMT, TrPot, RootFact, self.var.dz)
        
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

        # print self.var.th[0,0,0,0]
