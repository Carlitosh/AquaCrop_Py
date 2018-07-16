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
    
    def initial(self):

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
        self.var.th = th

    # def getState(self):
    #     result = {}
    #     state_vars_soil = [
    #         'AerDays','IrrCum','IrrNetCum',
    #         'DaySubmerged','Epot','Tpot','WTinSoil','AerDaysComp',
    #         'TrRatio','SurfaceStorage','SurfaceStorageIni','th',
    #         'Wsurf','EvapZ','Stage2','Wstage2']

    #     state_vars_crop = [
    #         # 'SeasonCounter',
    #         'DelayedGDDs','DelayedCDs','PctLagPhase','tEarlySen',
    #         'DAP','PreAdj','CropMature','CropDead','Germination',
    #         'PrematSenes','HarvestFlag','Stage','Fpre','Fpost',
    #         'fpost_dwn','fpost_upp','Fpol','sCor1','sCor2',
    #         'GrowthStage','CC','CCadj','CC_NS','CCadj_NS','B',
    #         'B_NS','HI','HIadj','CCxAct','CCxAct_NS','CCxW',
    #         'CCxW_NS','CCxEarlySen','rCor','Zroot','CC0adj']

    #     state_vars = state_vars_soil + state_vars_crop
    #     for var in state_vars:
    #         result[var] = vars(self.var)[var]

    def dynamic(self):
        # Condition to identify crops which are not being grown or crops which
        # have only just finished being grown. The water content of crops
        # meeting this condition is used to compute the area-weighted initial
        # condition
        if np.any(self.var.GrowingSeasonDayOne):
            cond1 = np.logical_not(self.var.GrowingSeasonIndex) | self.var.GrowingSeasonDayOne
            cond1 = np.broadcast_to(cond1, self.var.th.shape)
            th = np.copy(self.var.th)
            th[np.logical_not(cond1)] = np.nan
            th_ave = np.nanmean(th, axis=1)  # average along crop dimension
            th_ave = th_ave[:,None,:,:] * np.ones((self.var.nCrop))[None,:,None,None]
            cond2 = np.broadcast_to(self.var.GrowingSeasonDayOne, self.var.th.shape)
            self.var.th[cond2] = th_ave[cond2]
            
