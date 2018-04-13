#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
import shutil
import sys
import math
import gc

import pcraster as pcr

import virtualOS as vos
import meteo
import groundwater
import CropParameters as cropParams
import landSurface

import logging
logger = logging.getLogger(__name__)

class InitialConditions(object):

    def __init__(self, iniItems, currTimeStep, soilParams, cropParams, fieldMgmtParams):
        object.__init__(self)

        attr = vos.getMapAttributesALL(iniItems.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])
        self.nLayer = int(iniItems.soilOptions['nLayer'])
        self.nComp = int(iniItems.soilOptions['nComp'])
        
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))
        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))
                
        # Counters
        self.counter_names = ['AgeDays','AgeDays_NS','AerDays','IrrCum',
                              'IrrNetCum','DelayedGDDs','DelayedCDs',
                              'PctLagPhase','tEarlySen','GDDcum',
                              'DaySubmerged','DAP','Epot','Tpot']
        for nm in self.counter_names:
            vars(self)[nm] = arr_zeros
        
        # States
        self.state_names = ['WTinSoil','PreAdj','CropMature','CropDead',
                            'Germination','PrematSenes','HarvestFlag']
        for nm in self.state_names:
            vars(self)[nm] = arr_zeros.astype(bool)

        # Harvest index
        self.harvest_index_names_1 = ['Stage','Fpre','Fpost',
                                      'fpost_dwn','fpost_upp']
        for nm in self.harvest_index_names_1:
            vars(self)[nm] = arr_ones

        self.harvest_index_names_0 = ['HIcor_Asum','HIcor_Bsum','Fpol',
                                      'sCor1','sCor2']
        for nm in self.harvest_index_names_0:
            vars(self)[nm] = arr_zeros
            
        # Growth stage
        self.GrowthStage = arr_zeros

        # Aeration stress
        self.AerDaysComp = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))

        # Transpiration
        self.TrRatio = arr_ones
        
        # Crop growth
        self.crop_growth_names_0 = ['CC','CCadj','CC_NS','CCadj_NS','B','B_NS','HI','HIadj','CCxAct','CCxAct_NS','CCxW','CCxW_NS','CCxEarlySen','CCprev']
        for nm in self.crop_growth_names_zeros:
            vars(self)[nm] = arr_zeros

        self.rcor = arr_ones

        # If the start of the simulation equals the planting date of any crops
        # then set Zroot and CC0adj to Zmin and CC0 respectively. Otherwise set
        # to zero
        cond1 = (currTimeStep._startTime.timetuple().tm_yday == cropParams.PlantingDate)
        Zroot = np.zeros(cropParams.PlantingDate.shape)
        Zroot[cond1] = cropParams.Zmin[cond1]
        self.Zroot = np.max(Zroot, axis=0)

        CC0adj = np.zeros(cropParams.PlantingDate.shape)
        CC0adj[cond1] = cropParams.CC0[cond1]
        self.CC0adj = np.max(CC0adj, axis=0)

        # Surface storage between bunds
        cond2 = (fieldMgmtParams.Bunds == 0) & (fieldMgmtParams.zBund > 0.001)
        SurfaceStorage = np.zeros(fieldMgmtParams.Bunds.shape)
        SurfaceStorage[cond1] = fieldMgmtParams.BundWater[cond1]
        SurfaceStorage = np.clip(SurfaceStorage, None, fieldMgmtParams.zBund)
        self.SurfaceStorage = SurfaceStorage
        self.SurfaceStorageIni = SurfaceStorage
        
        # Define initial water contents
        self.initialConditionFileNC = iniItems.landSurfaceOptions['initialConditionNC']
        th = vos.netcdf2PCRobjCloneWithoutTime(
            self.initialConditionFileNC,'th', cloneMapFileName=self.cloneMap)
        init_cond_type = iniItems.globalOptions['initialConditionType']
        init_cond_interp_method = iniItems.globalOptions['initialConditionInterpMethod']
        
        if init_cond_type is 'Num':
            th = th
        elif init_cond_type is 'Pct':
            # NB original Matlab code allows users to supply initial conditions
            # at specified depth intervals or for specified layer. At the moment
            # we only allow the latter, although this may change in the future.
            th = self.soil.th_wp + ((th / 100) * (self.soil.th_fc - self.soil.th_wp))

        # NB original Matlab code also allows users to supply initial condition
        # as field capacity, wilting point or saturation. We do not explicitly
        # include this option but of course it is possible to do this by
        # supplying a netCDF file containing the values. However, note that in
        # the original version if field capacity is specified and groundwater
        # table is present the initial water content is set to th_fc_adj once
        # this variable has been computed.

        # Interpolate values to all soil compartments
        if init_cond_interp_method is 'Layer':
            self.th = th[self.soil.layerIndex,:]
        elif init_cond_interp_method is 'Depth':
            pass                # TODO

        # If groundwater table is present in soil profile then set all water
        # contents below the water table to saturation
        # Add dimensions to dz,dzsum
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = soilParams.dz[:,None,None,None] * arr_ones
        dzsum = soilParams.dzsum[:,None,None,None] * arr_ones

        # Water table in soil profile: calculate horizontal inflow; get
        # groundwater table elevation on current day
        zBot = np.cumsum(dz, axis=0)
        zTop = zBot - dz
        zMid = (zTop + zBot) / 2

        # For compartments below water table, set to saturation
        cond1 = (self.WTinSoil & (zMid >= soilParams.zGW) & (self.th < soilParams.th_s[layerIndex,:]))
        self.th[cond1] = soilParams.th_s[layerIndex,:][cond11]

    def reset(self, currTimeStep, cropParams):
        cond1 = (cropParams.CropSequence & (currTimeStep.doy == cropParams.PlantingDate))

        # Reset counters
        for nm in self.counter_names:
            vars(self)[nm][cond1] = 0

        # Reset states
        for nm in self.state_names:
            vars(self)[nm][cond1] = 0

        # Reset harvest index
        for nm in self.harvest_index_names_1:
            vars(self)[nm][cond1] = 1

        for nm in self.harvest_index_names_0:
            vars(self)[nm][cond1] = 0

        self.GrowthStage[cond1] = 0

        cond1_comp = np.broadcast_to(cond1, self.AerDaysComp.shape)
        self.AerDaysComp[cond1_comp] = 0
        
        self.TrRatio[cond1] = 1

        # Reset crop growth
        for nm in self.crop_growth_names_0:
            vars(self)[nm][cond1] = 0

        self.rCor[cond1] = 1
        self.CC0adj[cond1] = cropParams.CC0[cond1]
        self.Zroot[cond1] = cropParams.Zmin[cond1]

        # Update C02 concentration
        # NB this is done in CropParameters class
        
        # Reset soil water conditions (if not running off-season)
        # TODO

        # Update crop parameters (if in GDD mode)
        # TODO
        

        
            
        
            

        
