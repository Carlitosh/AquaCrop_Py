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
import scipy.interpolate as interpolate

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
        init_cond_interp_method = self.var._configuration.globalOptions['InterpMethod']
        # init_cond_type = self.var._configuration.globalOptions['initialConditionType']

        if init_cond_interp_method.lower() == 'depth':
            depths = vos.netcdfDim2NumPy(
                self.var.initialConditionFileNC,'depth')

            zBot = np.cumsum(self.var.dz)
            zTop = zBot - self.var.dz
            zMid = (zTop + zBot) / 2

            # add zero point
            if depths[0] > 0:
                depths = np.concatenate([[0.], depths])
                th = np.concatenate([th[0,...][None,:], th], axis=0)

            if depths[-1] < zBot[-1]:
                depths = np.concatenate([depths, [zBot[-1]]])
                th = np.concatenate([th, th[-1,...][None,:]], axis=0)

            # now interpolate to compartments
            f_thini = interpolate.interp1d(depths, th, axis=0, bounds_error=False, fill_value=(th[0,...], th[-1,...]))
            th = f_thini(zMid)
        else:
            th = th[self.var.layerIndex,:,:]
               
        # NB original Matlab code also allows users to supply initial condition
        # as field capacity, wilting point or saturation. We do not explicitly
        # include this option but of course it is possible to do this by
        # supplying a netCDF file containing the values. However, note that in
        # the original version if field capacity is specified and groundwater
        # table is present the initial water content is set to th_fc_adj once
        # this variable has been computed.

        th = np.broadcast_to(th, (self.var.nCrop, self.var.nComp, self.var.nLat, self.var.nLon))
        self.var.th = np.copy(th)        

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
            cond1 = np.broadcast_to(cond1[:,None,:,:], self.var.th.shape)
            th = np.copy(self.var.th)
            th[np.logical_not(cond1)] = np.nan
            th_ave = np.nanmean(th, axis=0)  # average along crop dimension
            th_ave = th_ave[None,:,:,:] * np.ones((self.var.nCrop))[:,None,None,None]
            cond2 = np.broadcast_to(self.var.GrowingSeasonDayOne[:,None,:,:], self.var.th.shape)
            self.var.th[cond2] = th_ave[cond2]
