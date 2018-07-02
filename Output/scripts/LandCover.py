#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
import shutil
import sys
import math
import gc
import numpy as np

from decimal import Decimal

import pcraster as pcr

from crop_growth_funs import *

import VirtualOS as vos
import Meteo as meteo
import CO2
import Groundwater as groundwater
import CropParameters as cropPars
import SoilAndTopoParameters as soilPars

import logging
logger = logging.getLogger(__name__)

class LandCover(object):

    def __init__(self, configuration, landmask, meteo, currTimeStep, initialState = None):

        self._configuration = configuration
        self._modelTime = currTimeStep

        # clone map, land mask
        pcr.setclone(configuration.cloneMap)
        self.cloneMap = configuration.cloneMap
        self.landmask = vos.readPCRmapClone(configuration.globalOptions['landmask'],
                                            configuration.cloneMap,
                                            configuration.tmpDir,
                                            configuration.globalOptions['inputDir'],
                                            True)
        attr = vos.getMapAttributesALL(self.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])

        # Load parameters
        # ###############
        self.crop_pars = cropPars.CropParameters(self._configuration, self.landmask)
        self.crop_pars.read()
        self.crop_pars.compute_variables(self._modelTime, meteo)

        self.nRotation = self.crop_pars.nRotation
        self.nCrop = self.crop_pars.nCrop
        
        # Set initial conditions
        # ######################
        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))
        self.DelayedGDDs = np.copy(arr_zeros)
        self.DelayedCDs  = np.copy(arr_zeros)
        self.PctLagPhase = np.copy(arr_zeros)
        self.tEarlySen   = np.copy(arr_zeros)
        self.GDDcum      = np.copy(arr_zeros)
        self.DAP         = np.copy(arr_zeros)

        self.AgeDays      = np.copy(arr_zeros)
        self.AgeDays_NS   = np.copy(arr_zeros)
        
        self.PreAdj      = np.copy(arr_zeros.astype(bool))
        self.CropMature  = np.copy(arr_zeros.astype(bool))
        self.CropDead    = np.copy(arr_zeros.astype(bool))
        self.Germination = np.copy(arr_zeros.astype(bool))
        self.PrematSenes = np.copy(arr_zeros.astype(bool))
        self.HarvestFlag = np.copy(arr_zeros.astype(bool))

        self.Stage = np.copy(arr_ones)
        self.Fpre = np.copy(arr_ones)
        self.Fpost = np.copy(arr_ones)
        self.fpost_dwn = np.copy(arr_ones)
        self.fpost_upp = np.copy(arr_ones)
        self.HIcor_Asum = np.copy(np.copy(arr_zeros))
        self.HIcor_Bsum = np.copy(np.copy(arr_zeros))
        self.Fpol = np.copy(np.copy(arr_zeros))
        self.sCor1 = np.copy(np.copy(arr_zeros))
        self.sCor2 = np.copy(np.copy(arr_zeros))
        
        self.GrowthStage = np.copy(np.copy(arr_zeros))

        self.CC = np.copy(np.copy(arr_zeros))
        self.CCadj = np.copy(np.copy(arr_zeros))
        self.CC_NS = np.copy(np.copy(arr_zeros))
        self.CCadj_NS = np.copy(np.copy(arr_zeros))
        self.B = np.copy(np.copy(arr_zeros))
        self.B_NS = np.copy(np.copy(arr_zeros))
        self.HI = np.copy(np.copy(arr_zeros))
        self.HIref = np.copy(np.copy(arr_zeros))
        self.HIadj = np.copy(np.copy(arr_zeros))
        self.CCxAct = np.copy(np.copy(arr_zeros))
        self.CCxAct_NS = np.copy(np.copy(arr_zeros))
        self.CCxW = np.copy(np.copy(arr_zeros))
        self.CCxW_NS = np.copy(np.copy(arr_zeros))
        self.CCxEarlySen = np.copy(np.copy(arr_zeros))
        self.CCprev = np.copy(np.copy(arr_zeros))        
        self.rCor = np.copy(arr_ones)
        self.Y = np.copy(np.copy(arr_zeros))
        
        # If the start of the simulation equals the planting date of any crops
        # then set Zroot and CC0adj to Zmin and CC0 respectively. Otherwise set
        # to zero
        # cond1 = (currTimeStep._startTime.timetuple().tm_yday == self.crop_pars.PlantingDate)
        # Zroot = np.zeros(self.crop_pars.PlantingDate.shape)
        # Zroot[cond1] = self.crop_pars.Zmin[cond1]
        # self.Zroot = np.max(Zroot, axis=0)

        # CC0adj = np.zeros(self.crop_pars.PlantingDate.shape)
        # CC0adj[cond1] = self.crop_pars.CC0[cond1]
        # self.CC0adj = np.max(CC0adj, axis=0)
        self.Zroot = np.copy(np.copy(arr_zeros))
        self.CC0adj = np.copy(np.copy(arr_zeros))
        
        # Declare other variables
        # #######################
        
        # Indexes
        self.SeasonCounter = np.zeros((self.nCrop, self.nLat, self.nLon))
        # self.GrowingSeason = np.zeros((self.nCrop, self.nRotation, self.nLat, self.nLon)).astype(bool)
        self.CropIndex = np.zeros((self.nRotation, self.nLat, self.nLon))
        self.GrowingSeasonIndex = np.zeros((self.nRotation, self.nLat, self.nLon)).astype(bool)
        
        # Harvest index
        self.YieldForm = np.copy(arr_zeros)
        self.HIref = np.copy(arr_zeros)
        self.HIt = np.copy(arr_zeros)
        self.PctLagPhase = np.copy(arr_zeros)

        # Crop growth
        self.GDD = np.copy(arr_zeros)
        self.Ksw_Exp = np.copy(arr_zeros)
        self.Ksw_Sto = np.copy(arr_zeros)
        self.Ksw_Sen = np.copy(arr_zeros)
        self.Ksw_Pol = np.copy(arr_zeros)
        self.Ksw_StoLin = np.copy(arr_zeros)
        self.Kst_Bio = np.copy(arr_zeros)
        self.Kst_PolH = np.copy(arr_zeros)
        self.Kst_PolC = np.copy(arr_zeros)

    def getState(self):
        result = {}
        state_vars = ['SeasonCounter',
                      'DelayedGDDs','DelayedCDs','PctLagPhase','tEarlySen',
                      'DAP','PreAdj','CropMature','CropDead','Germination',
                      'PrematSenes','HarvestFlag','Stage','Fpre','Fpost',
                      'fpost_dwn','fpost_upp','Fpol','sCor1','sCor2',
                      'GrowthStage','CC','CCadj','CC_NS','CCadj_NS','B',
                      'B_NS','HI','HIadj','CCxAct','CCxAct_NS','CCxW',
                      'CCxW_NS','CCxEarlySen','rCor','Zroot','CC0adj']
        for var in state_vars:
            result[var] = vars(self)[var]
        
    def getPseudoState(self):
        pass
        
    def update(self, meteo, CO2, currTimeStep):
        """Function to update parameters for current crop grown as well 
        as counters pertaining to crop growth
        """

        # Update crop parameters for currently grown crops
        # ################################################

        self.crop_pars.compute_water_productivity_adjustment_factor(CO2)
        self.crop_pars.update(currTimeStep, meteo)

        # Add a rotation dimension to GrowingSeason array and multiply it
        # by CropSequence, with the resulting array showing which crop (if any)
        # is currently grown in each rotation considered.
        GrowingSeason = self.crop_pars.GrowingSeason[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        GrowingSeason *= self.crop_pars.CropSequence  # crop,rotation,lat,lon
        
        # Lastly, check whether the crop has died or reached maturity, then
        # see which rotations have crops currently growing
        # GrowingSeason *= np.logical_not(self.CropDead | self.CropMature)
        self.GrowingSeasonIndex = np.any(GrowingSeason, axis=0)
        
        # Get index of crops currently grown
        CropIndex = (np.arange(0, self.nCrop)[:,None,None,None] * np.ones((self.nRotation, self.nLon, self.nLat))[None,:,:,:])
        CropIndex *= GrowingSeason
        self.CropIndex = np.max(CropIndex, axis=0).astype(int)
        
        # A CropIndex value of 0 currently means either that the crop is not
        # currently grown, or that the first crop (corresponding to index 0)
        # is grown. This confusion is handled in the next section by multiplying
        # the parameter value by GrowingSeason, such that it has a value of
        # zero if no crops are growing.

        # Select crop parameters for current day
        I,J,K = np.ogrid[:self.nRotation,:self.nLat,:self.nLon]
        for nm in self.crop_pars.crop_parameter_names:
            param = getattr(self.crop_pars, nm)
            param = param[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
            vars(self)[nm] = param[self.CropIndex,I,J,K] * self.GrowingSeasonIndex

        # GrowingSeasonDayOne is a logical array showing rotations for which
        # today is the start of a growing season
        self.GrowingSeasonDayOne = currTimeStep.doy == self.PlantingDate
        self.reset_initial_conditions(currTimeStep)

        self.GrowingSeasonIndex *= np.logical_not(self.CropDead | self.CropMature)
        
        # Update counters
        # ###############

        # Increment days after planting
        self.DAP[self.GrowingSeasonIndex] += 1
        self.DAP[np.logical_not(self.GrowingSeasonIndex)] = 0

        cond = (self.GrowingSeasonIndex & (self.DAP > self.MaxCanopyCD))
        self.AgeDays_NS[cond] = (self.DAP - self.MaxCanopyCD)[cond]
        
        # Increment growing degree days
        self.growing_degree_day(currTimeStep, meteo)
        self.GDDcum[self.GrowingSeasonIndex] += self.GDD[self.GrowingSeasonIndex]
        self.GDDcum[np.logical_not(self.GrowingSeasonIndex)] = 0

        # Adjust growth stage
        self.growth_stage()

    def reset_initial_conditions(self, currTimeStep):
        """Function to reset initial conditions at the start of a new
        growing season
        """
        cond = self.GrowingSeasonDayOne
        
        self.GrowthStage[cond] = 0
        self.DelayedGDDs[cond] = 0
        self.DelayedCDs[cond] = 0
        self.PctLagPhase[cond] = 0
        self.tEarlySen[cond] = 0
        self.GDDcum[cond] = 0
        self.DAP[cond] = 0

        self.AgeDays[cond] = 0
        self.AgeDays_NS[cond] = 0
        
        # States
        self.PreAdj[cond] = False
        self.CropMature[cond] = False
        self.CropDead[cond] = False
        self.Germination[cond] = False
        self.PrematSenes[cond] = False
        self.HarvestFlag[cond] = False

        # Harvest index
        self.Stage[cond] = 1
        self.Fpre[cond] = 1
        self.Fpost[cond] = 1
        self.fpost_dwn[cond] = 1
        self.fpost_upp[cond] = 1
        self.HIcor_Asum[cond] = 0
        self.HIcor_Bsum[cond] = 0
        self.Fpol[cond] = 0
        self.sCor1[cond] = 0
        self.sCor2[cond] = 0
        
        # Growth stage
        self.GrowthStage[cond] = 0

        # Reset crop growth
        self.CC[cond] = 0
        self.CCadj[cond] = 0
        self.CC_NS[cond] = 0
        self.CCadj_NS[cond] = 0
        self.B[cond] = 0
        self.B_NS[cond] = 0
        self.HI[cond] = 0
        self.HIadj[cond] = 0
        self.CCxAct[cond] = 0
        self.CCxAct_NS[cond] = 0
        self.CCxW[cond] = 0
        self.CCxW_NS[cond] = 0
        self.CCxEarlySen[cond] = 0
        self.CCprev[cond] = 0

        self.rCor[cond] = 1
        self.CC0adj[cond] = self.CC0[cond]
        self.Zroot[cond] = self.Zmin[cond]

    def growing_degree_day(self, currTimeStep, meteo):
        """Function to calculate number of growing degree days on 
        current day
        """
        tmax = meteo.tmax[None,:,:] * np.ones((self.nCrop))[:,None,None]  # TODO: should this be rotation?
        tmin = meteo.tmin[None,:,:] * np.ones((self.nCrop))[:,None,None]
        self.GDD = growing_degree_day(tmax, tmin, self.Tbase, self.Tupp, self.crop_pars.GDDmethod)
            
    def germination(self, soilwater):
        """Function to check if crop has germinated"""
        WcProp = water_content_affecting_germination(
            soilwater.th,
            soilwater.soil_pars.th_fc_comp,
            soilwater.soil_pars.th_wp_comp,
            soilwater.soil_pars.dz,
            soilwater.soil_pars.dzsum,
            soilwater.soil_pars.zGerm)
        
        # Check if water content is above germination threshold
        cond4 = (self.GrowingSeasonIndex & (WcProp >= self.GermThr) & (np.logical_not(self.Germination)))
        self.Germination[cond4] = True

        # Increment delayed growth time counters if germination is yet to occur
        cond5 = (self.GrowingSeasonIndex & (np.logical_not(self.Germination)))
        self.DelayedCDs[cond5] += 1
        self.DelayedGDDs[cond5] += self.GDD[cond5]

        # Update ageing days counter
        DAPadj = (self.DAP - self.DelayedCDs)
        cond6 = (DAPadj > self.MaxCanopyCD) & self.GrowingSeasonIndex
        self.AgeDays[cond6] = (DAPadj - self.MaxCanopyCD)[cond6]
        
        self.Germination[np.logical_not(self.GrowingSeasonIndex)] = False
        self.DelayedCDs[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.DelayedGDDs[np.logical_not(self.GrowingSeasonIndex)] = 0

    def growth_stage(self):
        """Function to calculate number of growing degree days on 
        current day
        """
        if self.crop_pars.CalendarType == 1:
            tAdj = self.DAP - self.DelayedCDs
        elif self.crop_pars.CalendarType == 2:
            tAdj = self.GDDcum - self.DelayedGDDs

        # Update growth stage
        cond1 = (self.GrowingSeasonIndex & (tAdj <= self.Canopy10Pct))
        cond2 = (self.GrowingSeasonIndex & np.logical_not(cond1) & (tAdj <= self.MaxCanopy))
        cond3 = (self.GrowingSeasonIndex & np.logical_not(cond1 | cond2) & (tAdj <= self.Senescence))
        cond4 = (self.GrowingSeasonIndex & np.logical_not(cond1 | cond2 | cond3) & (tAdj > self.Senescence))

        self.GrowthStage[cond1] = 1
        self.GrowthStage[cond2] = 2
        self.GrowthStage[cond3] = 3
        self.GrowthStage[cond4] = 4
        self.GrowthStage[np.logical_not(self.GrowingSeasonIndex)] = 0
        
    def root_development(self, groundwater, soilwater):
        """Function to calculate root zone expansion"""

        dZr, Zr = root_development(
            self.GrowingSeasonIndex,
            self.crop_pars.CalendarType,
            self.DAP, self.DelayedCDs, self.DelayedGDDs, self.GDD, self.GDDcum, self.Zmin, self.Zmax, self.PctZmin, self.Emergence, self.MaxRooting, self.fshape_r, self.fshape_ex, soilwater.TrRatio, self.Germination)
        
        # Get new rooting depth
        self.Zroot[self.GrowingSeasonIndex] = (self.Zroot + dZr)[self.GrowingSeasonIndex]

        # Adjust root depth if restrictive soil layer is present that limits
        # depth of root expansion
        cond11 = (self.GrowingSeasonIndex & (soilwater.soil_pars.zRes > 0))
        cond111 = (cond11 & (self.Zroot > soilwater.soil_pars.zRes))
        self.rCor[cond111] = np.divide(
            (2 * (self.Zroot / soilwater.soil_pars.zRes) * ((self.SxTop + self.SxBot) / 2) - self.SxTop),
            self.SxBot, out=np.zeros_like(Zr), where=cond111)[cond111]        
        self.Zroot[cond111] = soilwater.soil_pars.zRes[cond111]

        # Limit rooting depth if groundwater table is present (roots cannot
        # develop below the water table)
        if groundwater.WaterTable:
            zGW = np.copy(soilwater.zGW)
            cond12 = ((zGW > 0) & (self.Zroot > zGW))
            self.Zroot[cond12] = np.clip(zGW, self.Zmin, None)[cond12]

        # No root system outside of growing season
        self.Zroot[np.logical_not(self.GrowingSeasonIndex)] = 0

    def water_stress(self, meteo, soilwater, beta):
        """Function to calculate water stress coefficients"""        
        p_up = np.concatenate((self.p_up1[None,:], self.p_up2[None,:], self.p_up3[None,:], self.p_up4[None,:]), axis=0)
        p_lo = np.concatenate((self.p_lo1[None,:], self.p_lo2[None,:], self.p_lo3[None,:], self.p_lo4[None,:]), axis=0)
        fshape_w = np.concatenate((self.fshape_w1[None,:], self.fshape_w2[None,:], self.fshape_w3[None,:], self.fshape_w4[None,:]), axis=0)

        et0 = (meteo.referencePotET[None,:,:] * np.ones((self.nRotation))[:,None,None])
        self.Ksw_Exp, self.Ksw_Sto, self.Ksw_Sen, self.Ksw_Pol, self.Ksw_StoLin = water_stress(
            p_lo, p_up, fshape_w, et0, self.ETadj, self.tEarlySen,
            beta, self.beta, soilwater.TAW, soilwater.Dr)

    def canopy_cover(self, meteo, soilwater):
        """Function to simulate canopy growth/decline"""
        
        # Preallocate some variables
        CCxAdj = np.zeros((self.nRotation, self.nLat, self.nLon))
        CDCadj = np.zeros((self.nRotation, self.nLat, self.nLon))
        CCsen = np.zeros((self.nRotation, self.nLat, self.nLon))

        # Store initial condition
        self.CCprev = np.copy(self.CC)
        self.CC_NSprev = np.copy(self.CC_NS)
        
        # Calculate root zone water content, and determine if water stress is
        # occurring
        # self.root_zone_water()
        self.water_stress(meteo, soilwater, beta = True)

        # Get canopy cover growth over time
        if self.crop_pars.CalendarType == 1:
            tCC = np.copy(self.DAP)
            dtCC = np.ones((self.nRotation, self.nLat, self.nLon))
            tCCadj = self.DAP - self.DelayedCDs
        elif self.crop_pars.CalendarType == 2:
            tCC = np.copy(self.GDDcum)
            dtCC = np.copy(self.GDD)
            tCCadj = self.GDDcum - self.DelayedGDDs

        # Canopy development (potential)
        # ##############################
        self.CC_NS, self.CCxAct_NS, self.CCxW_NS = potential_canopy_development(
            self.GrowingSeasonIndex,
            self.Emergence, self.Maturity, self.Senescence,
            tCC, dtCC,
            self.CanopyDevEnd,
            self.CC0, self.CGC, self.CCx, self.CDC,
            self.CC_NS, self.CCxAct_NS, self.CCxW_NS)

        # Canopy development (actual)
        # ###########################
        self.CC, self.CC0adj, self.CCxAct = actual_canopy_development(
            self.GrowingSeasonIndex,
            self.Emergence, self.Maturity, self.Senescence,
            tCC, tCCadj, dtCC,
            self.CanopyDevEnd,
            self.CC0, self.CC0adj, self.CGC, self.CCx, self.CDC, self.Ksw_Exp, self.CC, self.CCprev, self.CCxAct)
        
        # Check for crop growth termination: if the following conditions are
        # met, the crop has died
        cond6 = (self.GrowingSeasonIndex & ((tCCadj <= self.Maturity) & (tCCadj >= self.Emergence) & (tCCadj > self.CanopyDevEnd)))
        cond63 = (cond6 & ((self.CC < 0.001) & np.logical_not(self.CropDead)))
        self.CC[cond63] = 0
        self.CropDead[cond63] = True

        # Canopy senescence due to water stress (actual)
        # ##############################################
        cond7 = (self.GrowingSeasonIndex & (tCCadj >= self.Emergence))

        # Check for early canopy senescence starting/continuing due to severe
        # water stress
        cond71 = (cond7 & ((tCCadj < self.Senescence) | (self.tEarlySen > 0)))

        # Early canopy senescence
        cond711 = (cond71 & (self.Ksw_Sen < 1))
        self.PrematSenes[cond711] = True

        # No prior early senescence
        cond7111 = (cond711 & (self.tEarlySen == 0))
        self.CCxEarlySen[cond7111] = self.CCprev[cond7111]

        # Increment early senescence GDD counter
        self.tEarlySen[cond711] += dtCC[cond711]

        update_CC_after_senescence(
            self.GrowingSeasonIndex, self.CC, self.CCprev, self.CC0, self.CC0adj, self.CGC, self.CDC, self.CCx, self.CCxAct, self.CCxW, self.CCxEarlySen, self.Ksw_Sen, self.tEarlySen, tCCadj, dtCC, self.Emergence, self.Senescence, self.PrematSenes, self.CropDead)

        # ##################################
        # adjust for micro-advective effects
        
        # Check to ensure potential CC is not slightly lower than actual
        cond8 = (self.GrowingSeasonIndex & (self.CC_NS < self.CC))
        self.CC_NS[cond8] = self.CC[cond8]

        cond81 = (cond8 & (tCC < self.CanopyDevEnd))
        self.CCxAct_NS[cond81] = self.CC_NS[cond81]

        # Actual (with water stress)
        self.CCadj[self.GrowingSeasonIndex] = ((1.72 * self.CC) - (self.CC ** 2) + (0.3 * (self.CC ** 3)))[self.GrowingSeasonIndex]

        # Potential (without water stress)
        self.CCadj_NS[self.GrowingSeasonIndex] = ((1.72 * self.CC_NS) - (self.CC_NS ** 2) + (0.3 * (self.CC_NS ** 3)))[self.GrowingSeasonIndex]

        # TODO: I don't think it is necessary to set these variables to zero here
        # No canopy outside growing season - set values to zero
        self.CC[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.CCadj[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.CC_NS[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.CCadj_NS[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.CCxW[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.CCxAct[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.CCxW_NS[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.CCxAct_NS[np.logical_not(self.GrowingSeasonIndex)] = 0

    def adjust_canopy_cover(self, soilwater):
        """Function to adjust canopy cover according to 
        transpiration
        """
        cond1 = (self.GrowingSeasonIndex & ((self.CC - self.CCprev) > 0.005) & (soilwater.TrAct > 0))
        self.CC[cond1] = self.CCprev[cond1]
        
    def harvest_index_ref_current_day(self):
        """Function to calculate reference (no adjustment for stress 
        effects) harvest index on current day
        """

        # Check if in yield formation period
        if self.crop_pars.CalendarType == 1:
            tAdj = self.DAP - self.DelayedCDs
        elif self.crop_pars.CalendarType == 2:
            tAdj = self.GDDcum - self.DelayedGDDs
        self.YieldForm = (self.GrowingSeasonIndex & (tAdj > self.HIstart))
        
        # Get time for harvest index calculation
        self.HIt = self.DAP - self.DelayedCDs - self.HIstartCD - 1

        harvest_index_ref_current_day(
            self.GrowingSeasonIndex, self.CropType,
            tAdj,
            self.HIt, self.HIref, self.HIini, self.HIGC, self.HI0, self.PctLagPhase,
            self.CCprev, self.CCmin, self.CCx, self.tLinSwitch, self.dHILinear)
        
    def temperature_stress(self, meteo):
        """Function to calculate temperature stress coefficients"""

        # Add rotation dimension to meteo vars
        tmin = meteo.tmin[None,:,:] * np.ones((self.nRotation))[:,None,None]
        tmax = meteo.tmax[None,:,:] * np.ones((self.nRotation))[:,None,None]
        
        # Calculate temperature stress coefficients affecting crop pollination
        self.Kst_Bio = temperature_stress_biomass(self.BioTempStress, self.GDD, self.GDD_up, self.GDD_lo)

        KsPol_up = 1
        KsPol_lo = 0.001
        self.Kst_PolH = temperature_stress_heat(tmax, self.Tmax_lo, self.Tmax_up, self.PolHeatStress, self.fshape_b, KsPol_up, KsPol_lo)
        self.Kst_PolC = temperature_stress_cold(tmin, self.Tmin_lo, self.Tmin_up, self.PolColdStress, self.fshape_b, KsPol_up, KsPol_lo)
        
    def biomass_accumulation(self, meteo, soilwater):
        """Function to calculate biomass accumulation"""

        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))
        
        et0 = meteo.referencePotET[None,:,:] * np.ones((self.nRotation))[:,None,None]
        self.temperature_stress(meteo)

        fswitch = np.zeros((self.nRotation, self.nLat, self.nLon))
        WPadj = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond1 = (self.GrowingSeasonIndex & (((self.CropType == 2) | (self.CropType == 3)) & (self.HIref > 0)))

        # # Adjust WP for reproductive stage
        # cond11 = (cond1 & (self.Determinant == 1))
        # fswitch[cond11] = (self.PctLagPhase / 100)[cond11]
        # cond12 = (cond1 & np.logical_not(cond11))
        # cond121 = (cond12 < (self.YldFormCD / 3))
        # fswitch[cond121] = np.divide(self.HIt, (self.YldFormCD / 3), out=np.copy(arr_zeros), where=self.YldFormCD!=0)[cond121]
        # cond122 = (cond12 & np.logical_not(cond121))
        # fswitch[cond122] = 1
        # WPadj[cond1] = (self.WP * (1 - (1 - self.WPy / 100) * fswitch))[cond1]
        # cond2 = (self.GrowingSeasonIndex & np.logical_not(cond1))
        # WPadj[cond2] = self.WP[cond2]
        WPadj = adjust_WP_for_reproductive_stage(
            self.GrowingSeasonIndex, self.CropType, self.HIref, self.HIt, self.PctLagPhase, self.Determinant, self.YldFormCD, self.WP, self.WPy)
        
        # Adjust WP for CO2 effects)
        WPadj *= self.fCO2

        # Calculate biomass accumulation on current day
        dB_NS = WPadj * (soilwater.TrPot_NS / et0) * self.Kst_Bio  # TODO: check correct TrPot is being used
        dB = WPadj * (soilwater.TrAct / et0) * self.Kst_Bio  # TODO: check correct TrAct is being used

        # Update biomass accumulation
        self.B += dB
        self.B_NS += dB_NS

        # No biomass accumulation outside growing season
        self.B[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.B_NS[np.logical_not(self.GrowingSeasonIndex)] = 0

    def harvest_index_adj_pre_anthesis(self):
        """Function to calculate adjustment to harvest index for 
        pre-anthesis water stress
        """
        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))

        cond0 = (self.GrowingSeasonIndex & self.YieldForm & (self.HIt >= 0) & ((self.CropType == 2) | (self.CropType == 3)) & np.logical_not(self.PreAdj))
        self.PreAdj[cond0] = True
        
        # Calculate adjustment
        Br = np.divide(self.B, self.B_NS, out=np.copy(arr_zeros), where=self.B_NS!=0)
        Br_range = np.log(self.dHI_pre, out=np.copy(arr_zeros), where=self.dHI_pre>0) / 5.62
        Br_upp = 1
        Br_low = 1 - Br_range
        Br_top = Br_upp - (Br_range / 3)

        # Get biomass ratio
        ratio_low_divd = (Br - Br_low)
        ratio_low_divs = (Br_top - Br_low)
        ratio_low = np.divide(ratio_low_divd, ratio_low_divs, out=np.copy(arr_zeros), where=ratio_low_divs!=0)
        ratio_upp_divd = (Br - Br_top)
        ratio_upp_divs = (Br_upp - Br_top)
        ratio_upp = np.divide(ratio_upp_divd, ratio_upp_divs, out=np.copy(arr_zeros), where=ratio_upp_divs!=0)

        # Calculate adjustment factor
        cond1 = (cond0 & ((Br >= Br_low) & (Br < Br_top)))
        self.Fpre[cond1] = (1 + (((1 + np.sin((1.5 - ratio_low) * np.pi)) / 2) * (self.dHI_pre / 100)))[cond1]
        cond2 = (cond0 & np.logical_not(cond1) & ((Br > Br_top) & (Br <= Br_upp)))
        self.Fpre[cond2] = (1 + (((1 + np.sin((0.5 + ratio_upp) * np.pi)) / 2) * (self.dHI_pre / 100)))[cond2]
        cond3 = (cond0 & np.logical_not(cond1 | cond2))
        self.Fpre[cond3] = 1

        # No green canopy left at start of flowering so no harvestable crop
        # will develop
        cond3 = (cond0 & (self.CC <= 0.01))
        self.Fpre[cond3] = 0
        
    def harvest_index_adj_pollination(self):
        """Function to calculate adjustment to harvest index for 
        failure of pollination due to water or temperature stress
        """
        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))
        FracFlow = np.copy(arr_zeros)
        t1 = np.copy(arr_zeros)
        t2 = np.copy(arr_zeros)
        F1 = np.copy(arr_zeros)
        F2 = np.copy(arr_zeros)
        F = np.copy(arr_zeros)

        cond0 = (self.GrowingSeasonIndex & self.YieldForm & (self.CropType == 3) & (self.HIt > 0) & (self.HIt <= self.FloweringCD))
        
        # Fractional flowering on previous day
        # cond1 = (HIt > 0)
        t1[cond0] = self.HIt[cond0] - 1
        cond11 = (cond0 & (t1 > 0))
        t1pct = 100 * np.divide(t1, self.FloweringCD, out=np.copy(arr_zeros), where=self.FloweringCD!=0)
        t1pct = np.clip(t1pct, 0, 100)
        F1[cond11] = (0.00558 * np.exp(0.63 * np.log(t1pct, out=np.copy(arr_zeros), where=t1pct>0)) - (0.000969 * t1pct) - 0.00383)[cond11]
        F1 = np.clip(F1, 0, None)

        # Fractional flowering on current day
        t2[cond0] = self.HIt[cond0]
        cond12 = (cond0 & (t2 > 0))
        t2pct = 100 * np.divide(t2, self.FloweringCD, out=np.copy(arr_zeros), where=self.FloweringCD!=0)
        t2pct = np.clip(t2pct, 0, 100)
        F2[cond12] = (0.00558 * np.exp(0.63 * np.log(t2pct, out=np.copy(arr_zeros), where=t2pct>0)) - (0.000969 * t2pct) - 0.00383)[cond12]
        F2 = np.clip(F2, 0, None)

        # Weight values
        cond13 = (cond0 & (np.abs(F1 - F2) >= 0.0000001))
        F[cond13] = (100 * np.divide(((F1 + F2) / 2), self.FloweringCD, out=np.copy(arr_zeros), where=self.FloweringCD!=0))[cond13]
        FracFlow[cond13] = F[cond13]
        
        # Calculate pollination adjustment for current day
        dFpol = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond2 = (cond0 & (self.CC >= self.CCmin))
        Ks = np.minimum(self.Ksw_Pol, self.Kst_PolC, self.Kst_PolH)
        dFpol[cond2] = (Ks * FracFlow * (1 + (self.exc / 100)))[cond2]

        # Calculate pollination adjustment to date
        self.Fpol += dFpol
        self.Fpol = np.clip(self.Fpol, None, 1)

    def harvest_index_adj_post_anthesis(self):
        """Function to calculate adjustment to harvest index for 
        post-anthesis water stress
        """
        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))

        cond0 = (self.GrowingSeasonIndex & self.YieldForm & (self.HIt > 0) & ((self.CropType == 2) | (self.CropType == 3)))
        
        # 1 Adjustment for leaf expansion
        tmax1 = self.CanopyDevEndCD - self.HIstartCD
        DAP = self.DAP - self.DelayedCDs
        cond1 = (cond0 & (DAP <= (self.CanopyDevEndCD + 1)) & (tmax1 > 0) & (self.Fpre > 0.99) & (self.CC > 0.001) & (self.a_HI > 0))
        dCor = (1 + np.divide((1 - self.Ksw_Exp), self.a_HI, out=np.copy(arr_zeros), where=self.a_HI!=0))
        self.sCor1[cond1] += np.divide(dCor, tmax1, out=np.copy(arr_zeros), where=tmax1!=0)[cond1]
        DayCor = (DAP - 1 - self.HIstartCD)
        self.fpost_upp[cond1] = (np.divide(tmax1, DayCor, out=np.copy(arr_zeros), where=DayCor!=0) * self.sCor1)[cond1]

        # 2 Adjustment for stomatal closure
        tmax2 = self.YldFormCD
        cond2 = (cond0 & (DAP <= (self.HIendCD + 1)) & (tmax2 > 0) & (self.Fpre > 0.99) & (self.CC > 0.001) & (self.b_HI > 0))
        dCor = ((np.exp(0.1 * np.log(self.Ksw_Sto))) * (1 - np.divide((1 - self.Ksw_Sto), self.b_HI, out=np.copy(arr_zeros), where=self.b_HI!=0)))
        self.sCor2[cond2] += np.divide(dCor, tmax2, out=np.copy(arr_zeros), where=tmax2!=0)[cond2]
        DayCor = (DAP - 1 - self.HIstartCD)
        self.fpost_dwn[cond2] = (np.divide(tmax2, DayCor, out=np.copy(arr_zeros), where=DayCor!=0) * self.sCor2)[cond2]

        # Determine total multiplier
        cond3 = (cond0 & (tmax1 == 0) & (tmax2 == 0))
        self.Fpost[cond3] = 1
        cond4 = (cond0 & np.logical_not(cond3))
        cond41 = (cond4 & (tmax2 == 0))
        self.Fpost[cond41] = self.fpost_upp[cond41]
        cond42 = (cond4 & (tmax1 <= tmax2) & np.logical_not(cond41))
        self.Fpost[cond42] = (self.fpost_dwn * np.divide(((tmax1 * self.fpost_upp) + (tmax2 - tmax1)), tmax2, out=np.copy(arr_zeros), where=tmax2!=0))[cond42]
        cond43 = (cond4 & np.logical_not(cond41 | cond42))
        self.Fpost[cond43] = (self.fpost_upp * np.divide(((tmax2 * self.fpost_dwn) + (tmax1 - tmax2)), tmax2, out=np.copy(arr_zeros), where=tmax2!=0))[cond43]

    def harvest_index(self, meteo, soilwater):
        """Function to simulate build up of harvest index"""
        
        # Calculate stresses, after updating root zone water content
        # self.root_zone_water()
        self.water_stress(meteo, soilwater, beta = True)

        # Calculate temperature stress
        self.temperature_stress(meteo)

        # Get reference harvest index on current day
        HIi = np.copy(self.HIref)

        # Calculate harvest index
        # #######################
        
        HIadj = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond1 = (self.GrowingSeasonIndex & self.YieldForm & (self.HIt >= 0))

        # Root/tuber or fruit/grain crops
        cond11 = (cond1 & ((self.CropType == 2) | (self.CropType == 3)))

        # # Determine adjustment for water stress before anthesis
        # cond111 = (cond11 & np.logical_not(self.PreAdj))
        # # self.PreAdj[cond111] = True
        self.harvest_index_adj_pre_anthesis()

        # Adjustment only for fruit/grain crops
        HImax = np.zeros((self.nRotation, self.nLat, self.nLon))  # TODO: is this in the right place?
        cond112 = (cond11 & (self.CropType == 3))
        # cond1121 = (cond112 & ((HIt > 0) & (HIt <= self.FloweringCD)))
        self.harvest_index_adj_pollination()
        HImax[cond112] = (self.Fpol * self.HI0)[cond112]
        cond113 = (cond11 & np.logical_not(cond112))
        HImax[cond113] = self.HI0[cond113]
                               
        # Determine adjustments for post-anthesis water stress
        # cond114 = (cond11 & (HIt > 0))
        self.harvest_index_adj_post_anthesis()

        # Limit HI to maximum allowable increase due to pre- and post-anthesis
        # water stress combinations
        HImult = np.ones((self.nRotation, self.nLat, self.nLon))
        HImult[cond11] = (self.Fpre * self.Fpost)[cond11]
        cond115 = (cond11 & (HImult > (1 + (self.dHI0 / 100))))
        HImult[cond115] = (1 + (self.dHI0 / 100))[cond115]

        # Determine harvest index on current day, adjusted for stress effects
        cond116 = (cond11 & (HImax >= HIi))
        HIadj[cond116] = (HImult * HIi)[cond116]
        
        cond117 = (cond11 & np.logical_not(cond116))
        HIadj[cond117] = (HImult * HImax)[cond117]

        # Leafy vegetable crops - no adjustment, harvest index equal to
        # reference value for current day
        cond12 = (cond1 & (self.CropType == 1))
        HIadj[cond12] = HIi[cond12]

        # Otherwise no build-up of harvest index if outside yield formation
        # period
        cond2 = (self.GrowingSeasonIndex & np.logical_not(cond1))
        HIi[cond2] = self.HI[cond2]
        HIadj[cond2] = self.HIadj[cond2]

        # Store final values for current time step
        self.HI[self.GrowingSeasonIndex] = HIi[self.GrowingSeasonIndex]
        self.HIadj[self.GrowingSeasonIndex] = HIadj[self.GrowingSeasonIndex]

        # No harvestable crop outside of a growing season
        self.HI[np.logical_not(self.GrowingSeasonIndex)] = 0
        self.HIadj[np.logical_not(self.GrowingSeasonIndex)] = 0

    def crop_yield(self):
        """Function to calculate crop yield"""
        cond1 = self.GrowingSeasonIndex
        self.Y[cond1] = ((self.B / 100) * self.HIadj)[cond1]
        cond11 = (cond1 & (((self.crop_pars.CalendarType == 1) & ((self.DAP - self.DelayedCDs) >= self.Maturity)) | ((self.crop_pars.CalendarType == 2) & ((self.GDDcum - self.DelayedGDDs) >= self.Maturity))))
        self.CropMature[cond11] = True
        self.Y[np.logical_not(self.GrowingSeasonIndex)] = 0
