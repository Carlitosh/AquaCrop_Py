#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
import shutil
import sys
import math
import gc
import numpy as np

import pcraster as pcr

from hydrology_funs import *

import VirtualOS as vos
import Meteo as meteo
import CO2
import Groundwater as groundwater
import CropParameters as cropParams
import SoilAndTopoParameters as soilParams
import IrrigationMgmtParameters as irriMgmtPars
import FieldMgmtParameters as fieldMgmtPars

import logging
logger = logging.getLogger(__name__)

class SoilWaterBalance(object):

    def __init__(self, configuration, landmask, groundwater, landcover, initialState = None):

        self._configuration = configuration
        # self._modelTime = currTimeStep

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
        self.nRotation = landcover.nRotation
        self.nCrop = landcover.nCrop
        
        # Load parameters
        # ###############        
        self.soil_pars = soilParams.SoilAndTopoParameters(self._configuration, self.landmask)
        self.soil_pars.read()
        self.nLayer = self.soil_pars.nLayer
        self.nComp = self.soil_pars.nComp
        
        self.irrigation_mgmt_pars = irriMgmtPars.IrrigationMgmtParameters(self._configuration, landmask)
        self.irrigation_mgmt_pars.read()

        self.field_mgmt_pars = fieldMgmtPars.FieldMgmtParameters(self._configuration, landmask)
        self.field_mgmt_pars.read()
                
        # Compute capillary rise parameters if water table is modelled
        # if groundwater.WaterTable:
        self.soil_pars.compute_capillary_rise_parameters()

        # Set initial conditions
        # ######################
        
        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))
        self.AerDays      = np.copy(arr_zeros)
        self.IrrCum       = np.copy(arr_zeros)
        self.IrrNetCum    = np.copy(arr_zeros)
        self.DaySubmerged = np.copy(arr_zeros)
        self.Epot         = np.copy(arr_zeros)
        self.Tpot         = np.copy(arr_zeros)
        
        self.WTinSoil     = np.copy(arr_zeros.astype(bool))

        # Aeration stress
        self.AerDaysComp  = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))

        # Transpiration
        self.TrRatio      = np.copy(arr_ones)
        
        # Surface storage between bunds
        cond1 = (self.field_mgmt_pars.Bunds == 0) & (self.field_mgmt_pars.zBund > 0.001)
        SurfaceStorage = np.zeros(self.field_mgmt_pars.Bunds.shape)
        SurfaceStorage[cond1] = self.field_mgmt_pars.BundWater[cond1]
        SurfaceStorage = np.clip(SurfaceStorage, None, self.field_mgmt_pars.zBund)
        self.SurfaceStorage = np.copy(SurfaceStorage)
        self.SurfaceStorageIni = np.copy(SurfaceStorage)

        # Define initial water contents
        self.initialConditionFileNC = self._configuration.globalOptions['initialConditionNC']
        th = vos.netcdf2PCRobjCloneWithoutTime(
            self.initialConditionFileNC,'th', cloneMapFileName=self.cloneMap)
        init_cond_type = self._configuration.globalOptions['initialConditionType']
        init_cond_interp_method = self._configuration.globalOptions['initialConditionInterpMethod']

        if init_cond_type is 'Num':
            th = th
        elif init_cond_type is 'Pct':
            # NB original Matlab code allows users to supply initial conditions
            # at specified depth intervals or for specified layer. At the moment
            # we only allow the latter, although this may change in the future.
            th = self.soil_pars.th_wp + ((th / 100) * (self.soil_pars.th_fc - self.soil_pars.th_wp))

        # NB original Matlab code also allows users to supply initial condition
        # as field capacity, wilting point or saturation. We do not explicitly
        # include this option but of course it is possible to do this by
        # supplying a netCDF file containing the values. However, note that in
        # the original version if field capacity is specified and groundwater
        # table is present the initial water content is set to th_fc_adj once
        # this variable has been computed.

        # # Interpolate values to all soil compartments
        # if init_cond_interp_method is 'Layer':
        #     self.th = th[self.soil_pars.layerIndex,:]
        # elif init_cond_interp_method is 'Depth':
        #     pass                # TODO

        # TODO: see above
        self.th = th            # TEMPORARY HACK

        # If groundwater table is present in soil profile then set all water
        # contents below the water table to saturation
        # Add dimensions to dz,dzsum
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = self.soil_pars.dz[:,None,None,None] * arr_ones
        dzsum = self.soil_pars.dzsum[:,None,None,None] * arr_ones

        # Water table in soil profile: calculate horizontal inflow; get
        # groundwater table elevation on current day
        zBot = np.cumsum(dz, axis=0)
        zTop = zBot - dz
        zMid = (zTop + zBot) / 2

        # For compartments below water table, set to saturation
        cond1 = (self.WTinSoil & (zMid >= groundwater.zGW) & (self.th < self.soil_pars.th_s[self.soil_pars.layerIndex,:]))
        self.th[cond1] = self.soil_pars.th_s[self.soil_pars.layerIndex,:][cond1]

        # Declare other variables
        # #######################

        # Groundwater
        self.zGW  = np.copy(arr_zeros)
        self.GwIn = np.copy(arr_zeros)
        self.th_fc_adj = self.soil_pars.th_fc[self.soil_pars.layerIndex,:]

        # Check for the presence of groundwater
        self.check_groundwater_table(groundwater)
        
        # Drainage
        self.FluxOut = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))
        self.DeepPerc = np.copy(arr_zeros)

        # Infiltration
        self.Runoff = np.copy(arr_zeros)
        self.Infl = np.copy(arr_zeros)

        # Capillary rise
        self.CrTot  = np.copy(arr_zeros)
        
        # Soil evaporation
        self.Wevap_Act = np.copy(arr_zeros)
        self.Wevap_Sat = np.copy(arr_zeros)
        self.Wevap_Fc  = np.copy(arr_zeros)
        self.Wevap_Wp  = np.copy(arr_zeros)
        self.Wevap_Dry = np.copy(arr_zeros)
        self.Stage2    = np.copy(arr_zeros.astype(bool))
        self.EvapZ     = np.copy(arr_zeros)
        self.Wstage2   = np.copy(arr_zeros)
        self.Wsurf     = np.copy(arr_zeros)
        
        # Irrigation
        self.Irr       = np.copy(arr_zeros)
        self.PreIrr    = np.copy(arr_zeros)
        self.IrrNet    = np.copy(arr_zeros)
        
        # Root zone
        self.thRZ_Act  = np.copy(arr_zeros)
        self.thRZ_Sat  = np.copy(arr_zeros)
        self.thRZ_Fc   = np.copy(arr_zeros)
        self.thRZ_Wp   = np.copy(arr_zeros)
        self.thRZ_Dry  = np.copy(arr_zeros)
        self.thRZ_Aer  = np.copy(arr_zeros)
        self.TAW       = np.copy(arr_zeros)
        self.Dr        = np.copy(arr_zeros)
        self.Wr        = np.copy(arr_zeros)

        # Aeration stress
        self.AerDaysComp = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))
                
        # Transpiration
        self.Ksa_Aer  = np.copy(arr_zeros)
        self.TrPot0   = np.copy(arr_zeros)
        self.TrPot_NS = np.copy(arr_zeros)
        self.TrAct    = np.copy(arr_zeros)
        self.TrAct0   = np.copy(arr_zeros)
        self.Tpot     = np.copy(arr_zeros)
        
    def getState(self):
        result = {}
        state_vars = ['AerDays','IrrCum','IrrNetCum',
                      'DaySubmerged','Epot','Tpot','WTinSoil','AerDaysComp',
                      'TrRatio','SurfaceStorage','SurfaceStorageIni','th',
                      'Wsurf','EvapZ','Stage2','Wstage2']
        # state_vars = ['AgeDays','AgeDays_NS','AerDays','IrrCum','IrrNetCum',
        #               'DaySubmerged','Epot','Tpot','WTinSoil','AerDaysComp',
        #               'TrRatio','SurfaceStorage','SurfaceStorageIni','th',
        #               'Wsurf','EvapZ','Stage2','Wstage2']
        for var in state_vars:
            result[var] = vars(self)[var]
        
    def getPseudoState(self):
        pass

    def update(self, landcover, currTimeStep):
        """Function to update irrigation and field management 
        parameters for currently grown crops
        """
        I,J,K = np.ogrid[:self.nRotation,:self.nLat,:self.nLon]
        for nm in self.irrigation_mgmt_pars.parameter_names:
            vars(self)[nm] = getattr(self.irrigation_mgmt_pars, nm)[landcover.CropIndex,I,J,K]

        for nm in self.field_mgmt_pars.parameter_names:
            param = getattr(self.field_mgmt_pars, nm)
            if nm in ['Bunds','zBund','BundWater']:
                vars(self)[nm] = param
            else:
                vars(self)[nm] = param[landcover.CropIndex,I,J,K]

        # reset initial conditions
        self.reset_initial_conditions(landcover, currTimeStep)

        cond1 = ((self.Bunds == 0) | (self.zBund < 0.001))
        self.DaySubmerged[cond1] = 0
                
    def reset_initial_conditions(self, landcover, currTimeStep):
        """Function to reset initial conditions at the start of a new
        growing season
        """    
        cond = landcover.GrowingSeasonDayOne
        self.AerDays[cond] = 0
        self.IrrCum[cond] = 0
        self.IrrNetCum[cond] = 0
        self.DaySubmerged[cond] = 0
        self.Epot[cond] = 0
        self.Tpot[cond] = 0
        
        self.WTinSoil[cond] = False
        cond_comp = np.broadcast_to(cond, self.AerDaysComp.shape)
        self.AerDaysComp[cond_comp] = 0        
        self.TrRatio[cond] = 1

    def check_groundwater_table(self, groundwater):
        """Function to check for presence of a groundwater table and, if 
        present, to adjust compartment water contents and field 
        capacities where necessary
        """
        if groundwater.WaterTable:

            # Copy depth to groundwater, and add rotation dimension for convenience
            self.zGW = groundwater.zGW[None,:,:] * np.ones((self.nRotation))[:,None,None]

            # get the mid point of each compartment
            zBot = np.cumsum(self.soil_pars.dz)
            zTop = zBot - self.soil_pars.dz
            zMid = (zTop + zBot) / 2
            zMid = zMid[:,None,None,None] * np.ones((self.nRotation,self.nLat,self.nLon))[None,:,:,:]

            # Check if water table is within modelled soil profile
            WTinSoilComp = (zMid >= self.zGW)
            self.th[WTinSoilComp] = self.soil_pars.th_s_comp[WTinSoilComp]

            # Flatten WTinSoilComp to provide an array with dimensions
            # (nrotation, nLat, nLon), indicating rotations where the water
            # table is in the soil profile
            self.WTinSoil = np.sum(WTinSoilComp, axis=0).astype(bool)

            # get Xmax
            Xmax = np.zeros((self.nComp,self.nRotation,self.nLat,self.nLon))
            cond1 = self.soil_pars.th_fc_comp <= 0.1
            cond2 = self.soil_pars.th_fc_comp >= 0.3
            cond3 = np.logical_not(cond1 | cond2) # i.e. 0.1 < fc < 0.3
            Xmax[cond1] = 1
            Xmax[cond2] = 2
            pF = 2 + 0.3 * (self.soil_pars.th_fc_comp - 0.1) / 0.2
            Xmax_cond3 = np.exp(pF * np.log(10)) / 100
            Xmax[cond3] = Xmax_cond3[cond3]
            cond4 = (self.zGW < 0) | ((self.zGW - zMid) >= Xmax)

            # Index of the compartment to which each element belongs (shallow ->
            # deep, i.e. 1 is the shallowest)
            compartment = (np.arange(1, self.nComp + 1)[:,None,None,None] * np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:])

            # Index of the lowest compartment (i.e. the maximum value) for which
            # cond4 is met, cast to all compartments (achieved by multiplying
            # compartments by cond4 to set elements that do not equal the
            # condition to zero, but retain the compartment number of elements
            # that do meet the condition
            cond4_max_compartment = (np.amax(compartment * cond4, axis=0)[None,:,:,:] * np.ones((self.nComp))[:,None,None,None])

            # Now, identify compartments that are shallower than the deepest
            # compartment for which cond4 is met
            cond4 = (compartment <= cond4_max_compartment)

            # 'cond4' is a special case because if ANY compartment meets the
            # condition then all overlying compartments are automatically assumed to
            # meet the condition. Thus in subsequent conditions we have to be careful
            # to ensure that True elements in 'cond4' do not also belong to 'cond5',
            # 'cond6' or 'cond7'. We use numpy.logical_not(...) for this purpose.
            cond5 = (self.soil_pars.th_fc_comp >= self.soil_pars.th_s_comp) & np.logical_not(cond4)
            cond6 = (zMid >= self.zGW) & np.logical_not(cond4 | cond5)
            cond7 = np.logical_not(cond4 | cond5 | cond6)
            dV = self.soil_pars.th_s_comp - self.soil_pars.th_fc_comp
            dFC = (dV / (Xmax ** 2)) * ((zMid - (self.zGW - Xmax)) ** 2)

            self.th_fc_adj[cond4] = self.soil_pars.th_fc_comp[cond4]
            self.th_fc_adj[cond5] = self.soil_pars.th_fc_comp[cond5]
            self.th_fc_adj[cond6] = self.soil_pars.th_s_comp[cond6]
            self.th_fc_adj[cond7] = self.soil_pars.th_fc_comp[cond7] + dFC[cond7]

        else:
            self.zGW = np.ones((self.nRotation, self.nLat, self.nLon)) * -999
            self.WTinSoil = np.full((self.nRotation, self.nLat, self.nLon), False)
            self.th_fc_adj = np.copy(self.soil_pars.th_fc_comp)

    def pre_irrigation(self, landcover):
        """Function to calculate pre-irrigation when in net irrigation 
        mode.
        """
        # Expand dz and dzsum to rotation, lat, lon
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = (self.soil_pars.dz[:,None,None,None] * arr_ones)
        dzsum = (self.soil_pars.dzsum[:,None,None,None] * arr_ones)

        # Calculate pre-irrigation requirement
        rootdepth = np.maximum(landcover.Zmin, landcover.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        thCrit = (self.soil_pars.th_wp_comp + ((self.NetIrrSMT / 100) * (self.soil_pars.th_fc_comp - self.soil_pars.th_wp_comp)))

        # Conditions for applying pre-irrigation
        cond1 = ((self.IrrMethod == 4) & (landcover.DAP == 1) & ((dzsum - dz) < rootdepth) & (self.th < thCrit))

        # Update pre-irrigation and root zone water content (mm)
        PreIrr_req = ((thCrit - self.th) * 1000 * dz)
        PreIrr_req[np.logical_not(cond1)] = 0
        self.PreIrr = np.sum(PreIrr_req, axis=0)

    def drainage(self):
        """Function to redistribute stored soil water"""        
        th, DeepPerc, FluxOut = drainage(self.th,
                                         self.soil_pars.th_s_comp,
                                         self.soil_pars.th_fc_comp,
                                         self.th_fc_adj,
                                         self.soil_pars.ksat_comp,
                                         self.soil_pars.tau_comp,
                                         self.soil_pars.dz,
                                         self.soil_pars.dzsum)
        self.th  = th
        self.DeepPerc = DeepPerc
        self.FluxOut = FluxOut
        
    def rainfall_partition(self, meteo):
        """Function to partition rainfall into surface runoff and 
        infiltration using the curve number approach.
        """
        Runoff, Infl = rainfall_partition(meteo.precipitation,
                                          self.th,
                                          self.soil_pars.th_fc_comp,
                                          self.soil_pars.th_wp_comp,
                                          self.soil_pars.zCN,
                                          self.soil_pars.AdjCN,
                                          self.soil_pars.CN,
                                          self.soil_pars.CNbot,
                                          self.soil_pars.CNtop,
                                          self.soil_pars.dz,
                                          self.soil_pars.dzsum,
                                          self.Bunds,
                                          self.zBund)
        self.Runoff = Runoff
        self.Infl = Infl
        
    def root_zone_water(self, landcover):
        """Function to calculate actual and total available water in the 
        root zone at current time step
        """
        thRZ, TAW, Dr, Wr = root_zone_water(self.th, self.soil_pars.th_s_comp, self.soil_pars.th_fc_comp, self.soil_pars.th_wp_comp, self.soil_pars.th_dry_comp, self.soil_pars.dz, landcover.Zmin, landcover.Zroot, landcover.Aer)
        self.thRZ_Act = thRZ['Act']
        self.thRZ_Sat = thRZ['Sat']
        self.thRZ_Fc  = thRZ['Fc']
        self.thRZ_Wp  = thRZ['Wp']
        self.thRZ_Dry = thRZ['Dry']
        self.thRZ_Aer = thRZ['Aer']
        self.TAW = TAW
        self.Dr = Dr
        self.Wr = Wr

    def irrigation(self, landcover, meteo):
        """Function to get irrigation depth for the current day"""

        SMT = np.concatenate((self.SMT1[None,:],
                              self.SMT2[None,:],
                              self.SMT3[None,:],
                              self.SMT4[None,:]), axis=0)

        # Determine adjustment for inflows and outflows on current day
        cond1 = (self.thRZ_Act > self.thRZ_Fc)
        rootdepth = np.maximum(landcover.Zmin, landcover.Zroot)
        AbvFc = ((self.thRZ_Act - self.thRZ_Fc) * 1000 * rootdepth)
        AbvFc[np.logical_not(cond1)] = 0
        WCadj = self.Tpot + self.Epot - meteo.precipitation + self.Runoff - AbvFc

        # Determine irrigation depth (mm/day) to be applied
        cond2 = (landcover.GrowingSeasonIndex & (self.IrrMethod == 0))
        self.Irr[cond2] = 0
        
        # If irrigation is based on soil moisture, get the soil moisture
        # target for the current growth stage and determine threshold to
        # initiate irrigation
        cond3 = (landcover.GrowingSeasonIndex & np.logical_not(cond2) & (self.IrrMethod == 1))
        I,J,K = np.ogrid[:self.nRotation,:self.nLat,:self.nLon]
        growth_stage_index = landcover.GrowthStage.astype(int) - 1
        SMT = SMT[growth_stage_index,I,J,K]
        IrrThr = np.round(((1 - SMT / 100) * self.TAW), 3)

        # If irrigation is based on a fixed interval, get number of days in
        # growing season so far (subtract 1 so that we always irrigate first
        # on day 1 of each growing season)
        cond4 = (landcover.GrowingSeasonIndex & (self.IrrMethod == 2))
        nDays = landcover.DAP - 1

        # Working on a copy, adjust depletion for inflows and outflows - same
        # for both soil moisture and interval based irrigation
        cond5 = (cond3 | cond4)
        Dr = np.copy(self.Dr)
        Dr[cond5] += WCadj[cond5]
        Dr[cond5] = np.clip(Dr, 0, None)[cond5]
        cond6 = ((cond3 & (Dr > IrrThr)) | (cond4 & ((nDays % self.IrrInterval) == 0)))        
        IrrReq = np.copy(Dr)
        EffAdj = ((100 - self.AppEff) + 100) / 100
        IrrReq *= EffAdj
        self.Irr[cond6] = np.clip(IrrReq, 0, self.MaxIrr)[cond6]
        cond7 = (cond5 & np.logical_not(cond6))
        self.Irr[cond7] = 0

        # If irrigation is based on a pre-defined schedule then the irrigation
        # requirement for each rotation is read from a netCDF file. Note that if
        # the option 'irrScheduleFileNC' is None, then nothing will be imported
        # and the irrigation requirement will be zero
        cond8 = (landcover.GrowingSeasonIndex & (self.IrrMethod == 3))
        if self.irrigation_mgmt_pars.irrScheduleFileNC != None:
            IrrReq = vos.netcdf2PCRobjClone(self.irrigation_mgmt_pars.irrScheduleFileNC,\
                                            "irrigation_schedule",\
                                            str(currTimeStep.fulldate),\
                                            useDoy = method_for_time_index,\
                                            cloneMapFileName = self.cloneMap,\
                                            LatitudeLongitude = True)
            self.Irr[cond8] = IrrReq[cond8]

        # Note that if irrigation is based on net irrigation then it is
        # performed after calculation of transpiration and hence is set to zero
        # at this point (not strictly necessary because Irr is initialized to
        # zero, but included for completeness)
        cond9 = (landcover.GrowingSeasonIndex & (self.IrrMethod == 4))
        self.Irr[cond9] = 0
        self.Irr[np.logical_not(landcover.GrowingSeasonIndex)] = 0
        self.IrrCum += self.Irr
        
    def infiltration(self):
        """Function to infiltrate incoming water (rainfall and 
        irrigation)
        """
        # Update infiltration rate for irrigation
        self.Infl += self.Irr * (self.AppEff / 100)

        SurfaceStorage, RunoffIni, ToStore = update_surface_storage(
            self.Bunds,
            self.zBund,
            self.soil_pars.ksat_comp,
            self.Infl,
            self.SurfaceStorage)
        
        Runoff, FluxOut, DeepPerc, thnew = infiltration(
            ToStore,
            self.FluxOut,
            self.th,
            self.soil_pars.th_s_comp,
            self.soil_pars.th_fc_comp,
            self.th_fc_adj,
            self.soil_pars.tau_comp,
            self.soil_pars.ksat_comp,
            self.soil_pars.dz)
        
        # Update total runoff
        Runoff += RunoffIni

        # Update surface storage (if bunds are present)
        self.SurfaceStorage = SurfaceStorage
        cond5 = ((Runoff > RunoffIni) & (self.Bunds == 1) & self.zBund > 0.001)
        self.SurfaceStorage[cond5] += (Runoff - RunoffIni)[cond5]

        # Limit surface storage to bund height: additional water above top of
        # bunds becomes runoff, and surface storage equals bund height
        cond51 = (cond5 & (self.SurfaceStorage > (self.zBund * 1000)))
        Runoff[cond51] = (RunoffIni + (self.SurfaceStorage - (self.zBund * 1000)))[cond51]
        self.SurfaceStorage[cond51] = (self.zBund * 1000)[cond51]
        cond52 = (cond5 & np.logical_not(cond51))
        Runoff[cond52] = RunoffIni[cond52]

        # Update water content, deep percolation, surface runoff, infiltration
        self.th = np.copy(thnew)
        self.DeepPerc += DeepPerc
        self.Infl -= Runoff
        self.Runoff += Runoff
        
    def capillary_rise(self, groundwater):
        """Function to calculate capillary rise from a shallow 
        groundwater table
        """

        if groundwater.WaterTable:

            # Get maximum capillary rise for bottom compartment
            zBot = np.sum(self.soil_pars.dz)              
            zBotMid = zBot - (self.soil_pars.dz[-1] / 2)  # depth to midpoint of bottom layer
            MaxCR = maximum_capillary_rise(
                self.soil_pars.ksat[-1,:], self.soil_pars.aCR[-1,:],
                self.soil_pars.bCR[-1,:], self.zGW, zBotMid)

            # Check for restrictions on upward flow caused by properties of
            # compartments that are not modelled in the soil water balance

            # Find top of next soil layer that is not within modelled soil profile
            zTopLayer = np.cumsum(self.soil_pars.zLayer[np.arange(0, (self.soil_pars.layerIndex[-1] + 1))])
            layeri = self.soil_pars.layerIndex[-1]  # layer number of bottom compartment
            LimCR = np.zeros((self.nRotation, self.nLat, self.nLon))

            # TODO: should this be nlayer, rather than nlayer-1
            while np.any(zTopLayer < self.zGW) & (layeri < (self.soil_pars.nLayer - 1)):
                layeri += 1
                LimCR = maximum_capillary_rise(
                    self.soil_pars.ksat[layeri,:], self.soil_pars.aCR[layeri,:],
                    self.soil_pars.bCR[layeri,:], self.zGW, zTopLayer)
                MaxCR = np.clip(MaxCR, None, LimCR)
                zTopLayer += self.soil_pars.zLayer[layeri]

            thnew, CrTot = capillary_rise_fun(
                self.th, self.soil_pars.th_fc_comp, self.th_fc_adj,
                self.soil_pars.th_wp_comp, self.soil_pars.fshape_cr_comp,
                self.soil_pars.aCR_comp, self.soil_pars.bCR_comp, MaxCR,
                self.FluxOut, self.zGW, self.soil_pars.dz)
            
            self.th = np.copy(thnew)
            self.CrTot = CrTot
                    
    def soil_evaporation(self, meteo, landcover, currTimeStep):
        """Function to calculate daily soil evaporation in AOS"""
        
        # Add rotation dimension to meteo vars
        et0 = meteo.referencePotET[None,:,:] * np.ones((self.nRotation))[:,None,None]
        prec = meteo.precipitation[None,:,:] * np.ones((self.nRotation))[:,None,None]

        # Prepare stage 2 evaporation (REW gone), if day one of simulation
        cond1 = (currTimeStep.timeStepPCR == 1)
        if currTimeStep.timeStepPCR == 1:
            self.Wsurf.fill(0)
            self.EvapZ = np.copy(self.soil_pars.EvapZmin)
            self.Wstage2 = prepare_stage_two_evaporation(
                self.th,
                self.soil_pars.th_s_comp,
                self.soil_pars.th_fc_comp,
                self.soil_pars.th_wp_comp,
                self.soil_pars.th_dry_comp,
                self.soil_pars.dz,
                self.soil_pars.REW,
                self.EvapZ)
            
        # Prepare soil evaporation stage 1 - adjust water in surface evaporation
        # layer for any infiltration (only do this if rainfall occurs or when
        # irrigation is triggered)
        cond2 = ((prec > 0) | np.any(((self.Irr > 0) & (self.IrrMethod != 4)), axis=0))
        cond21 = (cond2 & (self.Infl > 0))
        self.Wsurf[cond21] = self.Infl[cond21]
        self.Wsurf[cond21] = np.clip(self.Wsurf, None, self.soil_pars.REW)[cond21]
        self.Wstage2[cond21] = 0  # TODO: is this right?
        self.EvapZ[cond21] = self.soil_pars.EvapZmin[cond21]
        self.Stage2[cond21] = False

        # Calculate potential soil evaporation rate (mm/day)
        tAdj = adjust_time_for_delayed_development(
            landcover.crop_pars.CalendarType,
            landcover.DAP,
            landcover.DelayedCDs,
            landcover.DelayedGDDs,
            landcover.GDDcum)

        EsPot = potential_soil_evaporation_rate(
            landcover.GrowingSeasonIndex,
            tAdj,
            landcover.CC,
            landcover.CCadj,
            landcover.CCxAct,
            landcover.CCxW,
            landcover.Senescence,
            landcover.PrematSenes,
            et0,
            self.soil_pars.Kex,
            self.soil_pars.fwcc)
        
        # Adjust potential soil evaporation for mulches and/or partial wetting
        EsPotMul = adjust_potential_soil_evaporation_for_mulches(
            landcover.GrowingSeasonIndex,
            self.SurfaceStorage,
            EsPot,
            self.Mulches,
            self.fMulch,
            self.MulchPctGS,
            self.MulchPctOS)

        # Partial surface wetting by irrigation
        EsPotIrr = adjust_potential_soil_evaporation_for_irrigation(
            self.SurfaceStorage,
            EsPot,
            self.IrrMethod,
            self.Irr,
            prec,
            self.WetSurf)

        # Assign minimum value (mulches and partial wetting don't combine)
        EsPot = np.minimum(EsPotIrr, EsPotMul)

        # Initialise actual evaporation counter
        self.EsAct = np.zeros((self.nRotation, self.nLat, self.nLon))

        # Surface evaporation        
        EsActSurf = surface_evaporation(self.SurfaceStorage, EsPot)
        self.EsAct += EsActSurf
        cond = ((self.SurfaceStorage > 0) & (self.SurfaceStorage <= EsPot))
        self.SurfaceStorage -= EsActSurf
        self.Wsurf[cond] = self.soil_pars.REW[cond]
        self.Wstage2[cond] = 0
        self.EvapZ[cond] = self.soil_pars.EvapZmin[cond]
        
        # Stage 1 evaporation
        self.th, self.Wsurf, ToExtract, self.Wstage2 = soil_evaporation_stage_one(
            self.th,
            self.soil_pars.th_s_comp,
            self.soil_pars.th_fc_comp,
            self.soil_pars.th_wp_comp,
            self.soil_pars.th_dry_comp,
            EsPot,
            self.EsAct,
            self.soil_pars.dz,
            self.soil_pars.dzsum,
            self.soil_pars.EvapZmin,
            self.soil_pars.REW,
            self.EvapZ,
            self.Wsurf,
            self.Wstage2)

        # Stage 2 evaporation
        self.th, self.EsAct = soil_evaporation_stage_two(
            self.th,
            self.soil_pars.th_s_comp,
            self.soil_pars.th_fc_comp,
            self.soil_pars.th_wp_comp,
            self.soil_pars.th_dry_comp,
            EsPot,
            self.EsAct,
            self.soil_pars.dz,
            self.soil_pars.dzsum,
            self.soil_pars.EvapZmin,
            self.soil_pars.EvapZmax,
            self.soil_pars.REW,
            self.EvapZ,
            self.soil_pars.fWrelExp,
            self.soil_pars.fevap,
            self.Wstage2,
            ToExtract,
            EvapTimeSteps=20)
                
        # Store potential evaporation for irrigation calculations on next day
        self.Epot = np.copy(EsPot)

    def aeration_stress(self, landcover):
        """Function to calculate aeration stress coefficient"""

        # Determine aeration stress (root zone)
        cond1 = (self.thRZ_Act > self.thRZ_Aer)

        # Calculate aeration stress coefficient
        self.Ksa_Aer = np.ones((self.nRotation, self.nLat, self.nLon))
        cond11 = (cond1 & (self.AerDays < landcover.LagAer))
        x1 = (self.thRZ_Sat - self.thRZ_Act)
        x2 = (self.thRZ_Sat - self.thRZ_Aer)
        stress = (1 - np.divide(x1, x2, out=np.zeros_like(x1), where=x2!=0))
        self.Ksa_Aer[cond11] = (1 - ((self.AerDays / 3) * stress))[cond11]
        cond12 = (cond1 & np.logical_not(cond11))
        self.Ksa_Aer[cond12] = np.divide(x1, x2, out=np.zeros_like(x2), where=x2!=0)[cond12]

        # Increment aeration days counter, or set to zero if there is no stress
        self.AerDays[cond1] += 1
        self.AerDays[np.logical_not(cond1)] = 0
        self.AerDays = np.clip(self.AerDays, None, landcover.LagAer)

    def day_submerged(self, landcover):
        cond = (landcover.GrowingSeasonIndex & (self.SurfaceStorage > 0) & (self.DaySubmerged < landcover.LagAer))
        self.DaySubmerged[cond] += 1
        
    def transpiration(self, meteo, landcover, co2):
        """Function to calculate crop transpiration on current day"""

        # Add rotation dimension to ET0
        et0 = meteo.referencePotET[None,:,:] * np.ones((self.nRotation))[:,None,None]
        # arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))
        
        # potential transpiration
        # #######################
        
        # 1. No prior water stress
        self.TrPot_NS = potential_transpiration(
            landcover.GrowingSeasonIndex,
            landcover.Kcb,
            landcover.AgeDays_NS, landcover.fage, landcover.CCadj_NS, landcover.CCxW_NS, co2.conc, co2.RefConc, et0)
        
        # 2. Potential prior water stress and/or delayed development
        self.TrPot0 = potential_transpiration(
            landcover.GrowingSeasonIndex,
            landcover.Kcb,
            landcover.AgeDays, landcover.fage, landcover.CCadj, landcover.CCxW, co2.conc, co2.RefConc, et0)
        
        # Correct potential transpiration for dying green canopy effects
        cond9 = (landcover.GrowingSeasonIndex & (landcover.CC < landcover.CCxW))
        cond91 = (cond9 & (landcover.CCxW > 0.001) & (landcover.CC > 0.001))
        x = np.divide(landcover.CC, landcover.CCxW, out=np.zeros_like(landcover.CCxW), where=landcover.CCxW!=0)
        self.TrPot0[cond91] *= (x ** landcover.a_Tr)[cond91]

        # surface layer transpiration
        # ###########################
        
        # Potential transpiration counter
        TrPot = np.zeros((self.nRotation, self.nLat, self.nLon))
        TrPot_surf, self.TrAct0, self.AerDaysComp = surface_transpiration(
            landcover.GrowingSeasonIndex,
            self.SurfaceStorage,
            self.DaySubmerged,
            landcover.LagAer,
            self.AerDaysComp,
            self.TrPot0)
        TrPot += TrPot_surf

        # Update potential root zone transpiration for water stress
        # Maximum stress effect
        Ks = np.minimum(landcover.Ksw_StoLin, self.Ksa_Aer)
        
        # Update potential transpiration in root zone
        cond11 = (landcover.GrowingSeasonIndex & (self.IrrMethod != 4))
        TrPot[cond11] *= Ks[cond11]

        # Maximum sink term
        comp_sto, RootFact, SxComp = maximum_sink_term(
            landcover.GrowingSeasonIndex,
            self.IrrMethod,
            landcover.SxTop,
            landcover.SxBot,
            landcover.rCor,
            landcover.Zmin, landcover.Zroot, self.soil_pars.dz, self.soil_pars.dzsum)

        # Extract water
        self.th, self.AerDaysComp, self.TrAct = extract_transpiration_water(
            TrPot, landcover.GrowingSeasonIndex, comp_sto, self.DaySubmerged,
            self.IrrMethod, self.th, self.soil_pars.th_s_comp,
            self.soil_pars.th_fc_comp, self.soil_pars.th_wp_comp,
            self.soil_pars.th_dry_comp, landcover.ETadj, landcover.p_up2,
            landcover.p_lo2, landcover.fshape_w2, landcover.LagAer,
            landcover.Aer, self.AerDaysComp, et0, self.soil_pars.dz,
            SxComp, RootFact)
        
        # Add net irrigation water requirement
        # ####################################
        self.root_zone_water(landcover)
        self.th, self.IrrNet = add_net_irrigation(
            landcover.GrowingSeasonIndex,
            self.IrrMethod,
            self.th,
            self.soil_pars.th_wp_comp,
            self.soil_pars.th_fc_comp,
            self.thRZ_Wp, self.thRZ_Fc, self.thRZ_Act, self.NetIrrSMT, TrPot, RootFact, self.soil_pars.dz)
        
        # Update net irrigation counter for the growing season
        self.IrrNetCum += self.IrrNet
        
        # Add any surface transpiration to root zone total
        self.TrAct += self.TrAct0

        # Update transpiration ratio
        # ##########################
        
        cond17 = (landcover.GrowingSeasonIndex & (self.TrPot0 > 0))
        cond171 = (cond17 & (self.TrAct < self.TrPot0))
        self.TrRatio[cond171] = np.divide(self.TrAct, self.TrPot0, out=np.zeros_like(self.TrPot0), where=self.TrPot0!=0)[cond171]
        cond172 = (cond17 & np.logical_not(cond171))
        self.TrRatio[cond172] = 1
        cond18 = (landcover.GrowingSeasonIndex & np.logical_not(cond17))
        self.TrRatio[cond18] = 1
        self.TrRatio = np.clip(self.TrRatio, 0, 1)

        # No transpiration or irrigation if outside growing season
        # ########################################################
        
        self.TrAct[np.logical_not(landcover.GrowingSeasonIndex)] = 0
        self.TrPot0[np.logical_not(landcover.GrowingSeasonIndex)] = 0
        self.TrPot_NS[np.logical_not(landcover.GrowingSeasonIndex)] = 0
        self.IrrNet[np.logical_not(landcover.GrowingSeasonIndex)] = 0
        self.IrrNetCum[np.logical_not(landcover.GrowingSeasonIndex)] = 0

        # Store potential transpiration for irrigation calculations on next day
        self.Tpot = np.copy(self.TrPot0)

    def groundwater_inflow(self, groundwater):
        """Function to calculate capillary rise in the presence of a 
        shallow groundwater table
        """
        self.th, self.GwIn = groundwater_inflow(
            self.th, self.soil_pars.th_s_comp,
            self.WTinSoil,
            self.zGW,
            self.soil_pars.dz,
            self.soil_pars.dzsum)
        
    def add_pre_irrigation(self):
        self.IrrNet += self.PreIrr
