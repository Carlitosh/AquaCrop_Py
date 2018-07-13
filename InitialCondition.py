#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# AquaCrop
#
import os
import numpy as np
import pcraster as pcr
import VirtualOS as vos
import netCDF4 as nc
import datetime as datetime
import calendar as calendar

class InitialCondition(object):
    """Class to represent the initial condition of an AquaCrop run. 
    Although not yet implemented, this class should include a method 
    to read the model state from a dump netcdf, to enable a "warm" 
    start.
    """
    def __init__(self, InitialCondition_variable):
        self.var = InitialCondition_variable

class AquaCropInitialCondition(InitialCondition):
    
    def initial(self):

        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        arr_ones = np.ones((self.var.nRotation, self.var.nLat, self.var.nLon))

        # # Groundwater
        # self.var.WTinSoil = np.copy(arr_zeros.astype(bool))
        
        # # Irrigation
        # self.var.IrrCum = np.copy(arr_zeros)
        # self.var.IrrNetCum = np.copy(arr_zeros)

        # # Transpiration
        # self.var.AerDays = np.copy(arr_zeros)
        # self.var.AerDaysComp  = np.zeros((self.var.nComp, self.var.nRotation, self.var.nLat, self.var.nLon))
        # self.var.Tpot = np.copy(arr_zeros)        
        # self.var.TrRatio      = np.copy(arr_ones)
        # self.var.DaySubmerged = np.copy(arr_zeros)

        # # Soil evaporation
        # self.var.Epot = np.copy(arr_zeros)
        # self.var.Stage2 = np.copy(arr_zeros.astype(bool))
        # self.var.EvapZ = np.copy(arr_zeros)
        # self.var.Wstage2 = np.copy(arr_zeros)
        # self.var.Wsurf = np.copy(arr_zeros)

        # # Surface storage between bunds
        # cond1 = (self.var.Bunds == 0) & (self.var.zBund > 0.001)
        # SurfaceStorage = np.zeros(self.var.Bunds.shape)
        # SurfaceStorage[cond1] = self.var.BundWater[cond1]
        # SurfaceStorage = np.clip(SurfaceStorage, None, self.var.zBund)
        # self.var.SurfaceStorage = np.copy(SurfaceStorage)
        # self.var.SurfaceStorageIni = np.copy(SurfaceStorage)

        # Define initial water contents
        self.var.initialConditionFileNC = self.var._configuration.globalOptions['initialConditionNC']
        th = vos.netcdf2PCRobjCloneWithoutTime(
            self.var.initialConditionFileNC,'th', cloneMapFileName=self.var.cloneMap)
        init_cond_type = self.var._configuration.globalOptions['initialConditionType']
        init_cond_interp_method = self.var._configuration.globalOptions['initialConditionInterpMethod']

        if init_cond_type is 'Num':
            th = th
        elif init_cond_type is 'Pct':
            # NB original Matlab code allows users to supply initial conditions
            # at specified depth intervals or for specified layer. At the moment
            # we only allow the latter, although this may change in the future.
            th = self.var.th_wp + ((th / 100) * (self.var.th_fc - self.var.th_wp))

        # NB original Matlab code also allows users to supply initial condition
        # as field capacity, wilting point or saturation. We do not explicitly
        # include this option but of course it is possible to do this by
        # supplying a netCDF file containing the values. However, note that in
        # the original version if field capacity is specified and groundwater
        # table is present the initial water content is set to th_fc_adj once
        # this variable has been computed.

        # # Interpolate values to all soil compartments (TODO)
        # if init_cond_interp_method is 'Layer':
        #     self.th = th[self.soil_pars.layerIndex,:]
        # elif init_cond_interp_method is 'Depth':
        #     pass                # TODO
        self.var.th = th            # TEMPORARY HACK

        # Crop development
        # self.var.SeasonCounter = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))        
        # self.var.DelayedGDDs = np.copy(arr_zeros)
        # self.var.DelayedCDs  = np.copy(arr_zeros)
        # self.var.PctLagPhase = np.copy(arr_zeros)
        # self.var.tEarlySen   = np.copy(arr_zeros)
        self.var.GDDcum      = np.copy(arr_zeros)
        self.var.DAP         = np.copy(arr_zeros)
        # self.var.AgeDays      = np.copy(arr_zeros)
        # self.var.AgeDays_NS   = np.copy(arr_zeros)

        # Crop state
        self.var.PreAdj      = np.copy(arr_zeros.astype(bool))
        # self.var.CropMature  = np.copy(arr_zeros.astype(bool))
        # self.var.CropDead    = np.copy(arr_zeros.astype(bool))
        # self.var.Germination = np.copy(arr_zeros.astype(bool))
        # self.var.PrematSenes = np.copy(arr_zeros.astype(bool))
        self.var.HarvestFlag = np.copy(arr_zeros.astype(bool))  # NOT USED ANYWHERE

        self.var.Stage = np.copy(arr_ones)
        # self.var.Fpre = np.copy(arr_ones)
        # self.var.Fpost = np.copy(arr_ones)
        # self.var.fpost_dwn = np.copy(arr_ones)
        # self.var.fpost_upp = np.copy(arr_ones)
        # self.var.HIcor_Asum = np.copy(arr_zeros) # These vars do not appear to be used anywhere
        # self.var.HIcor_Bsum = np.copy(arr_zeros)
        # self.var.Fpol = np.copy(arr_zeros)
        # self.var.sCor1 = np.copy(arr_zeros)
        # self.var.sCor2 = np.copy(arr_zeros)
        
        # self.var.GrowthStage = np.copy(arr_zeros)
        # self.var.CC = np.copy(arr_zeros)
        # self.var.CCadj = np.copy(arr_zeros)
        # self.var.CC_NS = np.copy(arr_zeros)
        # self.var.CCadj_NS = np.copy(arr_zeros)
        # self.var.B = np.copy(arr_zeros)
        # self.var.B_NS = np.copy(arr_zeros)
        # self.var.HI = np.copy(arr_zeros)
        # self.var.HIref = np.copy(arr_zeros)  # TODO: is this used?
        # self.var.HIadj = np.copy(arr_zeros)
        # self.var.CCxAct = np.copy(arr_zeros)
        # self.var.CCxAct_NS = np.copy(arr_zeros)
        # self.var.CCxW = np.copy(arr_zeros)
        # self.var.CCxW_NS = np.copy(arr_zeros)
        # self.var.CCxEarlySen = np.copy(arr_zeros)
        # self.var.CCprev = np.copy(arr_zeros)  # TODO: is this used?
        # self.var.rCor = np.copy(arr_ones)
        # self.var.Y = np.copy(np.copy(arr_zeros))
        # self.var.Zroot = np.copy(arr_zeros)
        # self.var.CC0adj = np.copy(arr_zeros)
                
    def getState(self):
        result = {}
        state_vars_soil = [
            'AerDays','IrrCum','IrrNetCum',
            'DaySubmerged','Epot','Tpot','WTinSoil','AerDaysComp',
            'TrRatio','SurfaceStorage','SurfaceStorageIni','th',
            'Wsurf','EvapZ','Stage2','Wstage2']

        state_vars_crop = [
            # 'SeasonCounter',
            'DelayedGDDs','DelayedCDs','PctLagPhase','tEarlySen',
            'DAP','PreAdj','CropMature','CropDead','Germination',
            'PrematSenes','HarvestFlag','Stage','Fpre','Fpost',
            'fpost_dwn','fpost_upp','Fpol','sCor1','sCor2',
            'GrowthStage','CC','CCadj','CC_NS','CCadj_NS','B',
            'B_NS','HI','HIadj','CCxAct','CCxAct_NS','CCxW',
            'CCxW_NS','CCxEarlySen','rCor','Zroot','CC0adj']

        state_vars = state_vars_soil + state_vars_crop
        for var in state_vars:
            result[var] = vars(self.var)[var]
        
    # def update(self, landcover, currTimeStep):
    #     """Function to update irrigation and field management 
    #     parameters for currently grown crops
    #     """
    #     # reset initial conditions
    #     self.reset_initial_conditions(landcover, currTimeStep)

    #     cond1 = ((self.Bunds == 0) | (self.zBund < 0.001))
    #     self.DaySubmerged[cond1] = 0
                
    def reset_initial_conditions(self):
        """Function to reset initial conditions at the start of a new
        growing season
        """    
        cond = self.var.GrowingSeasonDayOne
        # self.var.AerDays[cond] = 0
        # self.var.IrrCum[cond] = 0
        # self.var.IrrNetCum[cond] = 0
        # self.var.DaySubmerged[cond] = 0
        # self.var.Epot[cond] = 0
        # self.var.Tpot[cond] = 0
        
        # # self.var.WTinSoil[cond] = False
        # cond_comp = np.broadcast_to(cond, self.var.AerDaysComp.shape)
        # self.var.AerDaysComp[cond_comp] = 0        
        # self.var.TrRatio[cond] = 1

        # Crop development
        # self.var.GrowthStage[cond] = 0
        # self.var.DelayedGDDs[cond] = 0
        # self.var.DelayedCDs[cond] = 0
        # self.var.PctLagPhase[cond] = 0
        # self.var.tEarlySen[cond] = 0
        self.var.GDDcum[cond] = 0
        self.var.DAP[cond] = 0
        # self.var.AgeDays[cond] = 0
        # self.var.AgeDays_NS[cond] = 0
        
        # Crop state
        # self.var.PreAdj[cond] = False
        # self.var.CropMature[cond] = False
        # self.var.CropDead[cond] = False
        # self.var.Germination[cond] = False
        # self.var.PrematSenes[cond] = False
        self.var.HarvestFlag[cond] = False  # NOT USED

        # Harvest index
        self.var.Stage[cond] = 1  # NOT USED
        # self.var.Fpre[cond] = 1
        # self.var.Fpost[cond] = 1
        # self.var.fpost_dwn[cond] = 1
        # self.var.fpost_upp[cond] = 1
        # self.var.HIcor_Asum[cond] = 0
        # self.var.HIcor_Bsum[cond] = 0
        # self.var.Fpol[cond] = 0
        # self.var.sCor1[cond] = 0
        # self.var.sCor2[cond] = 0

        # self.var.GrowthStage[cond] = 0  # TODO: does this make a difference?
        
        # Reset crop growth
        # self.var.CC[cond] = 0
        # self.var.CCadj[cond] = 0
        # self.var.CC_NS[cond] = 0
        # self.var.CCadj_NS[cond] = 0
        # self.var.B[cond] = 0
        # self.var.B_NS[cond] = 0
        # self.var.HI[cond] = 0
        # self.var.HIadj[cond] = 0
        # self.var.CCxAct[cond] = 0
        # self.var.CCxAct_NS[cond] = 0
        # self.var.CCxW[cond] = 0
        # self.var.CCxW_NS[cond] = 0
        # self.var.CCxEarlySen[cond] = 0
        # self.var.CCprev[cond] = 0

        # self.var.rCor[cond] = 1
        # self.var.CC0adj[cond] = self.var.CC0[cond]
        # self.var.Zroot[cond] = self.var.Zmin[cond]

    def growing_degree_day(self):
        tmax = self.var.tmax[None,:,:] * np.ones((self.var.nRotation))[:,None,None]
        tmin = self.var.tmin[None,:,:] * np.ones((self.var.nRotation))[:,None,None]
        if self.var.GDDmethod == 1:
            Tmean = ((tmax + tmin) / 2)
            Tmean = np.clip(Tmean, self.var.Tbase, self.var.Tupp)
        elif self.var.GDDmethod == 2:
            tmax = np.clip(tmax, self.var.Tbase, self.var.Tupp)
            tmin = np.clip(tmin, self.var.Tbase, self.var.Tupp)
            Tmean = ((tmax + tmin) / 2)
        elif self.var.GDDmethod == 3:
            tmax = np.clip(tmax, self.var.Tbase, self.var.Tupp)
            tmin = np.clip(tmin, None, self.var.Tupp)
            Tmean = np.clip(Tmean, self.var.Tbase, None)
        self.var.GDD = (Tmean - self.var.Tbase)
                
    def dynamic(self):

        # GrowingSeasonDayOne is a logical array showing rotations for which
        # today is the start of a growing season
        self.var.GrowingSeasonDayOne = self.var._modelTime.doy == self.var.PlantingDate
        self.reset_initial_conditions()
        self.var.GrowingSeasonIndex *= np.logical_not(self.var.CropDead | self.var.CropMature)

        # Necessary to do the following, because CropDead and CropMature aren't reset until later in the program running order (TODO)
        self.var.GrowingSeasonIndex[self.var.GrowingSeasonDayOne] = True

        # Update counters
        # Increment days after planting
        self.var.DAP[self.var.GrowingSeasonIndex] += 1
        self.var.DAP[np.logical_not(self.var.GrowingSeasonIndex)] = 0

        # cond = (self.var.GrowingSeasonIndex & (self.var.DAP > self.var.MaxCanopyCD))
        # self.var.AgeDays_NS[cond] = (self.var.DAP - self.var.MaxCanopyCD)[cond]
        
        # Increment growing degree days
        self.growing_degree_day()
        self.var.GDDcum[self.var.GrowingSeasonIndex] += self.var.GDD[self.var.GrowingSeasonIndex]
        self.var.GDDcum[np.logical_not(self.var.GrowingSeasonIndex)] = 0

        # Adjust growth stage
        self.var.growth_stage_module.dynamic()
                
        
class FAO56InitialCondition(InitialCondition):
    
    def initial(self):

        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        arr_ones = np.ones((self.var.nRotation, self.var.nLat, self.var.nLon))

        # Groundwater
        self.var.WTinSoil = np.copy(arr_zeros.astype(bool))
        
        # Irrigation
        self.var.IrrCum = np.copy(arr_zeros)
        self.var.IrrNetCum = np.copy(arr_zeros)

        # Transpiration/Soil evaporation (TODO)
        # self.var.AerDays = np.copy(arr_zeros)
        # self.var.AerDaysComp  = np.zeros((self.var.nComp, self.var.nRotation, self.var.nLat, self.var.nLon))
        # self.var.Tpot = np.copy(arr_zeros)        
        # self.var.TrRatio      = np.copy(arr_ones)
        self.var.ETpotCum = np.copy(arr_zeros)
        self.var.ETactCum = np.copy(arr_zeros)
        self.var.DaySubmerged = np.copy(arr_zeros)

        # # Soil evaporation
        # self.var.Epot = np.copy(arr_zeros)
        # self.var.Stage2 = np.copy(arr_zeros.astype(bool))
        # self.var.EvapZ = np.copy(arr_zeros)
        # self.var.Wstage2 = np.copy(arr_zeros)
        # self.var.Wsurf = np.copy(arr_zeros)

        # Surface storage between bunds
        cond1 = (self.var.Bunds == 0) & (self.var.zBund > 0.001)
        SurfaceStorage = np.zeros(self.var.Bunds.shape)
        SurfaceStorage[cond1] = self.var.BundWater[cond1]
        SurfaceStorage = np.clip(SurfaceStorage, None, self.var.zBund)
        self.var.SurfaceStorage = np.copy(SurfaceStorage)
        self.var.SurfaceStorageIni = np.copy(SurfaceStorage)

        # Define initial water contents
        self.var.initialConditionFileNC = self.var._configuration.globalOptions['initialConditionNC']
        th = vos.netcdf2PCRobjCloneWithoutTime(
            self.var.initialConditionFileNC,'th', cloneMapFileName=self.var.cloneMap)
        init_cond_type = self.var._configuration.globalOptions['initialConditionType']
        init_cond_interp_method = self.var._configuration.globalOptions['initialConditionInterpMethod']

        if init_cond_type is 'Num':
            th = th
        elif init_cond_type is 'Pct':
            # NB original Matlab code allows users to supply initial conditions
            # at specified depth intervals or for specified layer. At the moment
            # we only allow the latter, although this may change in the future.
            th = self.var.th_wp + ((th / 100) * (self.var.th_fc - self.var.th_wp))

        # NB original Matlab code also allows users to supply initial condition
        # as field capacity, wilting point or saturation. We do not explicitly
        # include this option but of course it is possible to do this by
        # supplying a netCDF file containing the values. However, note that in
        # the original version if field capacity is specified and groundwater
        # table is present the initial water content is set to th_fc_adj once
        # this variable has been computed.

        # # Interpolate values to all soil compartments (TODO)
        # if init_cond_interp_method is 'Layer':
        #     self.th = th[self.soil_pars.layerIndex,:]
        # elif init_cond_interp_method is 'Depth':
        #     pass                # TODO
        self.var.th = th            # TEMPORARY HACK

        # Crop development
        # self.var.SeasonCounter = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))        
        self.var.DAP         = np.copy(arr_zeros)
        # self.var.GrowthStage = np.copy(arr_zeros)
        self.var.Y = np.copy(np.copy(arr_zeros))
        self.var.Zroot = np.copy(arr_zeros)
                
    def getState(self):
        result = {}
        state_vars_soil = []
            # 'AerDays','IrrCum','IrrNetCum',
            # 'DaySubmerged','Epot','Tpot','WTinSoil','AerDaysComp',
            # 'TrRatio','SurfaceStorage','SurfaceStorageIni','th',
            # 'Wsurf','EvapZ','Stage2','Wstage2']

        state_vars_crop = []
            # # 'SeasonCounter',
            # 'DelayedGDDs','DelayedCDs','PctLagPhase','tEarlySen',
            # 'DAP','PreAdj','CropMature','CropDead','Germination',
            # 'PrematSenes','HarvestFlag','Stage','Fpre','Fpost',
            # 'fpost_dwn','fpost_upp','Fpol','sCor1','sCor2',
            # 'GrowthStage','CC','CCadj','CC_NS','CCadj_NS','B',
            # 'B_NS','HI','HIadj','CCxAct','CCxAct_NS','CCxW',
            # 'CCxW_NS','CCxEarlySen','rCor','Zroot','CC0adj']

        state_vars = state_vars_soil + state_vars_crop
        if len(state_vars) > 0:
            for var in state_vars:
                result[var] = vars(self.var)[var]
        
    # def update(self, landcover, currTimeStep):
    #     """Function to update irrigation and field management 
    #     parameters for currently grown crops
    #     """
    #     # reset initial conditions
    #     self.reset_initial_conditions(landcover, currTimeStep)

    #     cond1 = ((self.Bunds == 0) | (self.zBund < 0.001))
    #     self.DaySubmerged[cond1] = 0
                
    def reset_initial_conditions(self):
        """Function to reset initial conditions at the start of a new
        growing season
        """    
        cond = self.var.GrowingSeasonDayOne
        # self.var.AerDays[cond] = 0
        self.var.IrrCum[cond] = 0
        self.var.IrrNetCum[cond] = 0
        self.var.DaySubmerged[cond] = 0
        # self.var.Epot[cond] = 0
        # self.var.Tpot[cond] = 0
        self.var.ETpotCum[cond] = 0
        self.var.ETactCum[cond] = 0
        
        self.var.WTinSoil[cond] = False
        # cond_comp = np.broadcast_to(cond, self.var.AerDaysComp.shape)
        # self.var.AerDaysComp[cond_comp] = 0        
        # self.var.TrRatio[cond] = 1

        # # Groundwater
        # self.var.WTinSoil[cond] = False
        
        # # Irrigation
        # self.var.IrrCum[cond] = 0
        # self.var.IrrNetCum[cond] = 0

        # # Transpiration
        # self.var.AerDays[cond] = 0
        # cond_comp = np.broadcast_to(cond, self.var.AerDaysComp.shape)
        # self.var.AerDaysComp[cond_comp] = 0        
        # self.var.Tpot[cond] = 0
        # self.var.TrRatio[cond] = 1
        # self.var.DaySubmerged[cond] = 0

        # # Soil evaporation
        # self.var.Epot[cond] = 0 
        # self.var.Stage2[cond] = False
        # self.var.EvapZ[cond] = 0
        # self.var.Wstage2[cond] = 0
        # self.var.Wsurf[cond] = 0

        # Crop development
        # self.var.GrowthStage[cond] = 0
        # self.var.DelayedGDDs[cond] = 0
        # self.var.DelayedCDs[cond] = 0
        # self.var.PctLagPhase[cond] = 0
        # self.var.tEarlySen[cond] = 0
        # self.var.GDDcum[cond] = 0
        self.var.DAP[cond] = 0
        # self.var.AgeDays[cond] = 0
        # self.var.AgeDays_NS[cond] = 0
        
        # self.var.GrowthStage[cond] = 0  # TODO: does this make a difference?
        
        # Reset crop growth
        # self.var.Y[cond] = 0
        self.var.Zroot[cond] = self.var.Zmin[cond]

    def dynamic(self):

        # GrowingSeasonDayOne is a logical array showing rotations for which
        # today is the start of a growing season
        self.var.GrowingSeasonDayOne = self.var._modelTime.doy == self.var.PlantingDate
        self.reset_initial_conditions()
        # self.var.GrowingSeasonIndex *= np.logical_not(self.var.CropDead | self.var.CropMature)

        # Update counters
        # Increment days after planting
        self.var.DAP[self.var.GrowingSeasonIndex] += 1
        self.var.DAP[np.logical_not(self.var.GrowingSeasonIndex)] = 0

        # cond = (self.var.GrowingSeasonIndex & (self.var.DAP > self.var.MaxCanopyCD))
        # self.var.AgeDays_NS[cond] = (self.var.DAP - self.var.MaxCanopyCD)[cond]
        
        # # Increment growing degree days
        # self.growing_degree_day()
        # self.var.GDDcum[self.var.GrowingSeasonIndex] += self.var.GDD[self.var.GrowingSeasonIndex]
        # self.var.GDDcum[np.logical_not(self.var.GrowingSeasonIndex)] = 0

        # # Adjust growth stage
        # self.var.growth_stage_module.dynamic()
                
        
