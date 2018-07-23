#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class Infiltration(object):
    """Class to infiltrate incoming water (rainfall and irrigation)"""
    
    def __init__(self, Infiltration_variable):
        self.var = Infiltration_variable

    def initial(self):        
        # Surface storage between bunds
        cond1 = (self.var.Bunds == 0) & (self.var.zBund > 0.001)
        SurfaceStorage = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        SurfaceStorage[cond1] = self.var.BundWater[cond1]
        SurfaceStorage = np.clip(SurfaceStorage, None, self.var.zBund)
        self.var.SurfaceStorage = np.copy(SurfaceStorage)
        self.var.SurfaceStorageIni = np.copy(SurfaceStorage)

    def reset_initial_conditions(self):
        pass
        
    def dynamic(self):
        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()
            
        # dims = self.var.ksat_comp.shape
        # nc, nr, nlon, nlat = dims[0], dims[1], dims[2], dims[3]
        ToStore = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        RunoffIni = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        InflTot = self.var.Infl + self.var.SurfaceStorage
        thnew = np.copy(self.var.th)

        # Update infiltration rate for irrigation
        self.var.Infl += self.var.Irr * (self.var.AppEff / 100)

        # Infiltration limited by saturated hydraulic conductivity of surface
        # soil layer; additional water ponds on surface
        ksat_top = self.var.ksat_comp[:,0,...]
        cond1 = (self.var.Bunds == 1) & (self.var.zBund > 0.001)
        cond11 = (cond1 & (InflTot > 0))
        cond111 = (cond11 & (InflTot > ksat_top))
        ToStore[cond111] = ksat_top[cond111]
        self.var.SurfaceStorage[cond111] = (InflTot - ksat_top)[cond111]

        # Otherwise all water infiltrates and surface storage becomes zero
        cond112 = (cond11 & np.logical_not(cond111))
        ToStore[cond112] = InflTot[cond112]
        self.var.SurfaceStorage[cond112] = 0

        # Calculate additional RunoffIni if water overtops bunds
        cond113 = (cond11 & (self.var.SurfaceStorage > (self.var.zBund * 1000)))
        RunoffIni[cond113] = (self.var.SurfaceStorage - (self.var.zBund * 1000))[cond113]
        self.var.SurfaceStorage[cond113] = (self.var.zBund * 1000)[cond113]

        # Otherwise excess water does not overtop bunds and there is no RunoffIni
        cond114 = (cond11 & np.logical_not(cond113))
        RunoffIni[cond114] = 0

        # If total infiltration is zero then there is no storage or RunoffIni
        cond12 = (cond1 & np.logical_not(cond11))
        ToStore[cond12] = 0
        RunoffIni[cond12] = 0

        # If there are no bunds then infiltration is divided between RunoffIni
        # and infiltration according to saturated conductivity of surface
        # layer
        cond2 = (self.var.Bunds == 0)
        cond21 = (cond2 & (self.var.Infl > ksat_top))
        ToStore[cond21] = ksat_top[cond21]
        RunoffIni[cond21] = (self.var.Infl - ksat_top[0,:])[cond21]
        cond22 = (cond2 & np.logical_not(cond21))
        ToStore[cond22] = self.var.Infl[cond22]
        RunoffIni[cond22] = 0

        # #######################################################################

        # Initialize counters
        comp = 0
        Runoff = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        cond3_ini = (ToStore > 0)

        while (np.any(ToStore > 0) & (comp < self.var.nComp)):

            # Update condition
            cond3 = (ToStore > 0)

            # Calculate saturated drainage ability, drainage factor and
            # required drainage ability
            dthdtS = self.var.tau_comp[:,comp,...] * (self.var.th_s_comp[:,comp,...] - self.var.th_fc_comp[:,comp,...])
            factor = self.var.ksat_comp[:,comp,...] / (dthdtS * 1000 * self.var.dz[comp])
            dthdt0 = ToStore / (1000 * self.var.dz[comp])

            # Initialize water content for current layer
            theta0 = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

            # Check drainage ability
            cond31 = (cond3 & (dthdt0 < dthdtS))

            # Calculate water content needed to meet drainage dthdt0
            cond311 = (cond31 & (dthdt0 <= 0))
            theta0[cond311] = self.var.th_fc_adj[:,comp,...][cond311]
            cond312 = (cond31 & np.logical_not(cond311))
            A = (1 + ((dthdt0 * (np.exp(self.var.th_s_comp[:,comp,...] - self.var.th_fc_comp[:,comp,...]) - 1)) / (self.var.tau_comp[:,comp,...] * (self.var.th_s_comp[:,comp,...] - self.var.th_fc_comp[:,comp,...]))))
            theta0[cond312] = (self.var.th_fc_comp[:,comp,...] + np.log(A))[cond312]

            # Limit thX to between saturation and field capacity
            cond313 = (cond31 & (theta0 > self.var.th_s_comp[:,comp,...]))
            theta0[cond313] = self.var.th_s_comp[:,comp,...][cond313]
            cond314 = (cond31 & np.logical_not(cond313) & (theta0 < self.var.th_fc_adj[:,comp,...]))
            theta0[cond314] = self.var.th_fc_adj[:,comp,...][cond314]
            dthdt0[cond314] = 0

            # Limit water content and drainage to saturation
            cond32 = (cond3 & np.logical_not(cond31))
            theta0[cond32] = self.var.th_s_comp[:,comp,...][cond32]
            dthdt0[cond32] = dthdtS[cond32]

            # Calculate maximum water flow through compartment and total
            # drainage
            drainmax = factor * dthdt0 * 1000 * self.var.dz[comp]
            drainage = drainmax + self.var.FluxOut[:,comp,...]

            # Limit drainage to saturated hydraulic conductivity
            cond33 = (cond3 & (drainage > self.var.ksat_comp[:,comp,...]))
            drainmax[cond33] = (self.var.ksat_comp[:,comp,...] - self.var.FluxOut[:,comp,...])[cond33]

            # Line 117 of AOS_Infiltration.m
            # Calculate difference between threshold and current water contents
            diff = theta0 - self.var.th[:,comp,...]

            cond34 = (cond3 & (diff > 0))
            thnew[:,comp,...][cond34] += (ToStore / (1000 * self.var.dz[comp]))[cond34]

            # Remaining water that can infiltrate to compartments below
            cond341 = (cond34 & (thnew[:,comp,...] > theta0))
            ToStore[cond341] = ((thnew[:,comp,...] - theta0) * 1000 * self.var.dz[comp])[cond341]
            thnew[:,comp,...][cond341] = theta0[cond341]

            # Otherwise all infiltrating water has been stored
            cond342 = (cond34 & np.logical_not(cond341))
            ToStore[cond342] = 0

            # Update outflow from current compartment (drainage + infiltration
            # flows)
            self.var.FluxOut[:,comp,...][cond3] += ToStore[cond3]

            # Calculate backup of water into compartments above, and update
            # water to store
            excess = np.clip((ToStore - drainmax), 0, None)
            ToStore -= excess

            # Redistribute excess to compartments above
            precomp = comp + 1
            while (np.any(cond3 & (excess > 0)) & (precomp != 0)):

                # Update condition
                cond35 = (cond3 & (excess > 0))

                # Keep storing in compartments above until soil surface is
                # reached
                precomp -= 1

                # Update outflow from compartment
                self.var.FluxOut[:,precomp,...][cond35] = (self.var.FluxOut[:,precomp,...] - excess)[cond35]

                # Update water content and limit to saturation
                thnew[:,precomp,...][cond35] += (excess / (1000 * self.var.dz[precomp]))[cond35]
                cond351 = (cond35 & (thnew[:,precomp,...] > self.var.th_s_comp[:,precomp,...]))
                excess[cond351] = ((thnew[:,precomp,...] - self.var.th_s_comp[:,precomp,...]) * 1000 * self.var.dz[precomp])[cond351]
                thnew[:,precomp,...][cond351] = self.var.th_s_comp[:,precomp,...][cond351]
                cond352 = (cond35 & np.logical_not(cond351))
                excess[cond352] = 0

            # Any leftover water not stored becomes Runoff
            cond36 = (cond3 & (excess > 0))
            Runoff[cond36] += excess[cond36]

            # update comp
            comp += 1

        # Infiltration left to store after bottom compartment becomes deep
        # percolation (mm)
        DeepPerc = np.copy(ToStore)

        # Otherwise if ToStore equals zero there is no infiltration
        cond4 = np.logical_not(cond3_ini)
        DeepPerc[cond4] = 0
        Runoff[cond4] = 0
        
        # #######################################################################

        # Update total runoff
        Runoff += RunoffIni

        # Update surface storage (if bunds are present)
        # self.var.SurfaceStorage = SurfaceStorage
        cond5 = ((Runoff > RunoffIni) & (self.var.Bunds == 1) & self.var.zBund > 0.001)
        self.var.SurfaceStorage[cond5] += (Runoff - RunoffIni)[cond5]

        # Limit surface storage to bund height: additional water above top of
        # bunds becomes runoff, and surface storage equals bund height
        cond51 = (cond5 & (self.var.SurfaceStorage > (self.var.zBund * 1000)))
        Runoff[cond51] = (RunoffIni + (self.var.SurfaceStorage - (self.var.zBund * 1000)))[cond51]
        self.var.SurfaceStorage[cond51] = (self.var.zBund * 1000)[cond51]
        cond52 = (cond5 & np.logical_not(cond51))
        Runoff[cond52] = RunoffIni[cond52]

        # Update water content, deep percolation, surface runoff, infiltration
        self.var.th = np.copy(thnew)
        self.var.DeepPerc += DeepPerc
        self.var.Infl -= Runoff
        self.var.Runoff += Runoff        
