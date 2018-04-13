#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
import shutil
import sys
import math
import gc

import pcraster as pcr

import VirtualOS as vos
import Meteo as meteo
import C02
import Groundwater as groundwater
import CropParameters as cropParams
import SoilAndTopoParameters as soilParams
import InitialConditions as initCond

import logging
logger = logging.getLogger(__name__)

class AquaCrop(object):

    def __init__(self, configuration, currTimeStep, initialState = None):

        self._configuration = configuration
        self._modelTime = currTimeStep

        # clone map, land mask
        pcr.setclone(configuration.cloneMap)
        self.landmask = vos.readPCRmapClone(configuration.globalOptions['landmask'],\
                                            configuration.cloneMap,\
                                            configuration.tmpDir,\
                                            configuration.globalOptions['inputDir'],\
                                            True)

        # Soil parameters
        self.soil = soilParams.SoilAndTopoParameters(iniItems, self.landmask)
        self.soil.read()

        # Crop parameters
        self.crop = cropParams.CropParameters(iniItems, landmask)
        self.crop.read()
        self.crop.compute_variables()

        # Irrigation management parameters
        self.irrigation_mgmt = irriMgmtParams.IrrigationMgmtParameters(self._configuration, landmask)
        self.irrigation_mgmt.read()

        # Field management parameters
        self.field_mgmt = fieldMgmtParams.FieldMgmtParameters(self._configuration, landmask)
        self.field_mgmt.read()
        
        # Prepare sub-models
        self.meteo = meteo.Meteo(
            self._configuration,
            self.landmask,
            initialState)

        self.C02 = CO2.C02(
            self._configuration,
            self.landmask,
            initialState)
        
        self.groundwater = groundwater.Groundwater(
            self._configuration,
            self.landmask,
            initialState)

        # Compute capillary rise parameters if water table is modelled
        if groundwater.WaterTable:
            self.soil.compute_capillary_rise_parameters()

        # If groundwater is present, this function calculates adjusted FC
        self.soil.check_groundwater_table(groundwater)
        
        # Define variables
        # TODO: take out of function
        self.set_initial_conditions()

    @property
    def configuration(self):
        return self._configuration

    def set_initial_conditions(self):
        self.init_conditions = initCond.InitialConditions(self._configuration, self._modelTime, self.crop)
            
        # Groundwater
        self.GwIn = arr_zeros
        self.zGW = arr_zeros
        self.th_fc_adj = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))

        # Drainage
        self.FluxOut = arr_zeros
        self.DeepPerc = arr_zeros

        # Infiltration
        self.Runoff = arr_zeros
        self.Infl = arr_zeros
        self.SurfaceStorage = arr_zeros

        # Soil evaporation
        self.Wevap = {'Act': arr_zeros,
                      'Sat': arr_zeros,
                      'Fc': arr_zeros,
                      'Wp': arr_zeros,
                      'Dry': arr_zeros}
        self.Stage2 = arr_zeros.astype(bool)
        self.EvapZ = arr_zeros
        self.Wstage2 = arr_zeros
        self.Wsurf = arr_zeros
        
        # Irrigation
        self.Irr
        self.PreIrr_req
        self.PreIrr
        self.IrrNet = arr_zeros
        
        # Root zone
        self.thRZ = {'Act': arr_zeros,
                     'Sat': arr_zeros,
                     'Fc': arr_zeros,
                     'Wp': arr_zeros,
                     'Dry': arr_zeros,
                     'Aer': arr_zeros}
        self.TAW = arr_zeros
        self.Dr = arr_zeros
        
        # Harvest index
        self.YieldForm = arr_zeros
        self.HIref = arr_zeros
        self.PctLagPhase = arr_zeros

        # Aeration stress
        self.AerDaysComp = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))

        # Transpiration
        self.Ksa_Aer = arr_zeros
        self.TrPot0 = arr_zeros
        self.TrPot_NS = arr_zeros
        self.TrAct,
        self.TrAct0 = arr_zeros
        self.Tpot = arr_zeros
        
        # Crop growth
        self.GDD = arr_zeros
        self.Ksw = {'Exp': arr_zeros,
                    'Sto': arr_zeros,
                    'Sen': arr_zeros,
                    'Pol': arr_zeros,
                    'StoLin': arr_zeros}
        self.Kst = {'Bio': arr_zeros,
                    'PolH': arr_zeros,
                    'PolC': arr_zeros}

    # def reset_initial_conditions(self):
    #     pass
        
    def dumpState(self, outputDirectory, specific_date_string = None):
        pass

    def resume(self):
        pass

    def setState(self):
        logger.error("cannot set state")

    def getState(self):
        pass

    def getPseudoState(self):
        pass

    def getAllState(self):
        pass

    # def read_forcings(self):
    #     logger.info("Reading forcings for time %s", self._modelTime)
    #     self.meteo.read_forcings(self._modelTime)
    #     self.groundwater.read_forcings(self._modelTime)

    def update(self, meteo, currTimeStep):
        """Function to update parameters for current crop grown as well 
        as counters pertaining to crop growth
        """

        # Read forcings for current time step
        logger.info("Updating model for time %s", self._modelTime)
        self.meteo.read_forcings(self._modelTime)
        self.groundwater.read_forcings(self._modelTime)
        self.C02.read_forcings(self._modelTime)

        # Update season counter for crops planted on current day
        cond1 = (self.CropSequence & (currTimeStep.doy == self.PlantingDate))
        self.initCond.reset()   # TODO
        self.SeasonCounter[cond1] += 1

        # Update growing season
        cond2 = (self.SeasonCounter > 0)
        cond3 = ((self.PlantingDate <= currTimeStep.doy) & (currTimeStep.doy <= self.HarvestDate))
        cond4 = (PlantingDate > HarvestDate)
        cond5 = ((PlantingDate <= currTimeStep.doy) | (currTimeStep.doy <= HarvestDate))
        self.GrowingSeason = ((cond2 & cond3) | (cond2 & cond4 & cond5))

        # Get index of crops currently grown
        I,J,K = np.ogrid[:self.nRotation,:self.nLat,:self.nLon]
        crop_index = (
            np.arange(0, self.nCrop)[:,None,None,None]
            * np.ones((self.nRotation, self.nLon, self.nLat))[None,:,:,:])
        crop_index *= self.GrowingSeason
        crop_index = np.max(crop_index, axis=0).astype(int)
        growing_season_index = np.any(self.GrowingSeason, axis=0)

        # Increment days after planting
        self.DAP[growing_season_index] += 1
        self.DAP[np.logical_not(growing_season_index)] = 0

        # Increment growing degree days
        self.growing_degree_day(meteo)
        self.GDDcum[growing_season_index] += self.GDD[growing_season_index]
        self.GDDcum[np.logical_not(growing_season_index)] = 0

        # Select parameters for current day
        for nm in self.crop.parameter_names:
            param = getattr(self.crop, nm)
            param = param[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
            vars(self)[nm] = param[crop_index,I,J,K]

        for nm in self.irrigation_mgmt.parameter_names:
            vars(self)[nm] = getattr(self.irrigation_mgmt, nm)[crop_index,I,J,K]

        for nm in self.field_mgmt.parameter_names:
            vars(self)[nm] = getattr(self.field_mgmt, nm)[crop_index,I,J,K]
        
    # def check_groundwater_table(self, groundwater, currTimeStep):
    #     """Function to check for presence of a groundwater table and, if 
    #     present, to adjust compartment water contents and field 
    #     capacities where necessary
    #     """
    #     # expand soil properties to compartments
    #     th_fc = self.soil.th_fc[self.soil.layerIndex,:]
    #     th_s  = self.soil.th_s[self.soil.layerIndex,:]
        
    #     # These are options which apply to the entire study area (not just a
    #     # selection of cells), so we can use an if statement.
    #     if groundwater.WaterTable:
        
    #         # get the mid point of each compartment (TODO: put this in initialization)
    #         zBot = np.cumsum(self.soil.dz)
    #         zTop = zBot - self.soil.dz
    #         zMid = (zTop + zBot) / 2
        
    #         # Convenient to add compartment + rotation dimensions to groundwater
    #         # and rotation, latitude and longitude dimensions to zMid
    #         self.zGW = groundwater.zGW[None,:,:] * np.ones((self.nRotation))[:,:,None,None]
    #         zMid = zMid[:,None,None,None] * np.ones((self.nRotation,self.nLat,self.nLon))[None,:,:,:]

    #         # Check if water table is within modelled soil profile
    #         WTinSoilComp = (zMid >= self.zGW)
    #         self.th[WTinSoilComp] = th_s[WTinSoilComp]

    #         # Flatten WTinSoilComp to provide an array with dimensions
    #         # (nrotation, nLat, nLon), indicating rotations where the water
    #         # table is in the soil profile
    #         self.WTinSoil = np.sum(WTinSoilComp, axis=0).astype(bool)

    #         # get Xmax (TODO: find out what this variable represents)
    #         Xmax = np.zeros((self.soil.nComp,self.nRotation,self.nLat,self.nLon))
    #         cond1 = th_fc <= 0.1
    #         cond2 = th_fc >= 0.3
    #         cond3 = np.logical_not(cond1 | cond2) # i.e. 0.1 < fc < 0.3
    #         Xmax[cond1] = 1
    #         Xmax[cond2] = 2
    #         pF = 2 + 0.3 * (th_fc - 0.1) / 0.2
    #         Xmax_cond3 = np.exp(pF * np.log(10)) / 100
    #         Xmax[cond3] = Xmax_cond3[cond3]
    #         cond4 = (self.zGW < 0) | ((self.zGW - zMid) >= Xmax)

    #         # Index of the compartment to which each element belongs (shallow ->
    #         # deep, i.e. 1 is the shallowest)
    #         compartment = (
    #             np.arange(1, self.nComp + 1)[:,None,None,None]
    #             * np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:])

    #         # Index of the lowest compartment (i.e. the maximum value) for which
    #         # cond4 is met, cast to all compartments (achieved by multiplying
    #         # compartments by cond4 to set elements that do not equal the
    #         # condition to zero, but retain the compartment number of elements
    #         # that do meet the condition
    #         cond4_max_compartment = (
    #             np.amax(compartment * cond4, axis=0)[None,:,:,:]
    #             * np.ones((self.soil.nComp))[:,None,None,None])

    #         # Now, identify compartments that are shallower than the deepest
    #         # compartment for which cond4 is met
    #         cond4 = (compartment <= cond4_max_compartment)

    #         # 'cond4' is a special case because if ANY compartment meets the
    #         # condition then all overlying compartments are automatically assumed to
    #         # meet the condition. Thus in subsequent conditions we have to be careful
    #         # to ensure that True elements in 'cond4' do not also belong to 'cond5',
    #         # 'cond6' or 'cond7'. We use numpy.logical_not(...) for this purpose.
    #         cond5 = (th_fc >= th_s) & np.logical_not(cond4)
    #         cond6 = (zMid >= self.zGW) & np.logical_not(cond4 | cond5)
    #         cond7 = np.logical_not(cond4 | cond5 | cond6)
    #         dV = th_s - th_fc
    #         dFC = (dV / (Xmax ** 2)) * ((zMid - (self.zGW - Xmax)) ** 2)

    #         self.th_fc_adj[cond4] = th_fc[cond4]
    #         self.th_fc_adj[cond5] = th_fc[cond5]
    #         self.th_fc_adj[cond6] = th_s[cond6]
    #         self.th_fc_adj[cond7] = th_fc[cond7] + dFC[cond7]

    #     else:
    #         self.zGW = np.ones((self.nRotation, self.nLat, self.nLon)) * -999
    #         self.WTinSoil = np.full((self.nRotation, self.nLat, self.nLon), False)
    #         self.th_fc_adj = th_fc
    
    def pre_irrigation(self, rotation):
        """Function to calculate pre-irrigation when in net irrigation 
        mode.
        """
        # Expand soil properties to compartments
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]

        # Expand dz and dzsum to rotation, lat, lon
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = (self.soil.dz[:,None,None,None] * arr_ones)
        dzsum = (self.soil.dzsum[:,None,None,None] * arr_ones)

        # Calculate pre-irrigation requirement
        rootdepth = np.maximum(self.Zmin, self.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        thCrit = (th_wp + ((self.NetIrrSMT / 100) * (th_fc - th_wp)))

        # Conditions for applying pre-irrigation
        cond1 = ((self.IrrMethod == 4)
                 & (self.DAP == 1)
                 & ((dzsum - dz) < rootdepth)
                 & (self.th < thCrit))

        # Update pre-irrigation and root zone water content (mm)
        PreIrr_req = ((thCrit - self.th) * 1000 * dz)
        PreIrr_req[np.logical_not(cond1)] = 0
        self.PreIrr_req = PreIrr_req
        self.PreIrr = np.sum(PreIrr_req, axis=0)

    def drainage(self):
        """Function to redistribute stored soil water"""
        
        # Expand soil properties to compartments
        ksat = self.soil.ksat[self.soil.layerIndex,:]
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_s = self.soil.th_s[self.soil.layerIndex,:]
        tau = self.soil.tau[self.soil.layerIndex,:]
        
        # Preallocate arrays        
        drainsum = np.zeros((self.nRotation, self.nLat, self.nLon))
        thnew = self.th

        for comp in range(self.soil.nComp):

            # Calculate drainage ability of compartment ii
            dthdt = np.zeros((self.nRotation, self.nLat, self.nLon))
            cond1 = (self.th[comp,:] <= self.soil.th_fc_adj[comp,:])
            dthdt[cond1] = 0

            cond2 = (np.logical_not(cond1) & (self.th[comp,:] >= th_s[comp,:]))
            dthdt[cond2] = ((tau[comp,:] * th_s[comp,:]) - th_fc[comp,:])[cond2]

            cond3 = np.logical_not(cond1 | cond2)
            dthdt[cond3] = (
                tau[comp,:]
                * (th_s[comp,:] - th_fc[comp,:])
                * ((np.exp(self.th[comp,:] - th_fc[comp,:]) - 1)
                   / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond3]

            cond4 = (
                (cond2 | cond3)
                & ((self.th[comp,:] - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond4] = (self.th[comp,:] - self.soil.th_fc_adj[comp,:])[cond4]

            # Drainage from compartment ii (mm)
            draincomp = dthdt * self.soil.dz[comp]

            # Check drainage ability of compartment ii against cumulative
            # drainage from compartments above
            excess = np.zeros((self.nRotation, self.nLat, self.nLon))
            prethick = self.soil.dzsum[comp] - self.soil.dz[comp]
            drainmax = dthdt * 1000 * prethick
            drainability = (drainsum <= drainmax)

            # Drain compartment
            cond5 = drainability
            thnew[comp,:][cond5] = (self.th[comp,:] - dthdt)[cond5]
            
            # Update cumulative drainage (mm), restrict to saturated hydraulic
            # conductivity and adjust excess drainage flow
            drainsum[cond5] += draincomp[cond5]
            cond51 = (cond5 & (drainsum > ksat[comp,:]))
            excess[cond51] += (drainsum - ksat[comp,:])[cond51]
            drainsum[cond51] = ksat[comp,:][cond51]

            # ##################################################################
            # Calculate value of theta (thX) needed to provide a drainage
            # ability equal to cumulative drainage

            cond6 = np.logical_not(drainability)
            dthdt[cond6] = (drainsum / (1000 * prethick))[cond6]

            thX = np.zeros((self.nRotation, self.nLat, self.nLon))

            cond61 = (cond6 & (dthdt <= 0))
            thX[cond61] = self.soil.th_fc_adj[comp,:][cond61]
            
            cond62 = (cond6 & np.logical_not(cond61) & (tau[comp,:] > 0))
            A = (1 + ((dthdt * (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1))
                      / (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]))))
            thX[cond62] = (self.soil.th_fc_adj[comp,:] + np.log(A))[cond62]

            cond621 = (cond62 & (thX < self.soil.th_fc_adj[comp,:]))
            thX[cond621] = self.soil.th_fc_adj[comp,:][cond621]

            cond63 = (cond6 & np.logical_not(cond61 | cond62))
            thX[cond63] = (th_s[comp,:] + 0.01)[cond63]

            # ##################################################################
            # Check thX against hydraulic properties of current soil layer

            # Increase compartment ii water content with cumulative drainage
            cond64 = (np.logical_not(drainability) & (thX <= th_s[comp,:]))
            thnew[comp,:][cond64] = (
                self.th[comp,:]
                + (drainsum / (1000 * self.soil.dz[comp])))[cond64]

            # Check updated water content against thX
            cond641 = (cond64 & (thnew[comp,:] > thX))

            # Cumulative drainage is the drainage difference between theta_x and
            # new theta plus drainage ability at theta_x
            drainsum[cond641] = (
                (thnew[comp,:] - thX)
                * 1000 * self.soil.dz[comp])[cond641]
            
            # Calculate drainage ability for thX
            cond6411 = (cond641 & (thX <= self.soil.th_fc_adj[comp,:]))
            dthdt[cond6411] = 0
            cond6412 = (
                cond641
                & np.logical_not(cond6411)
                & (thX >= th_s[comp,:]))
            dthdt[cond6412] = (
                tau[comp,:]
                * (th_s[comp,:] - th_fc[comp,:]))[cond6412]
            cond6413 = (cond641 & np.logical_not(cond6411 | cond6412))
            dthdt[cond6413] = (
                tau[comp,:]
                * (th_s[comp,:] - th_fc[comp,:])
                * ((np.exp(thX - th_fc[comp,:]) - 1)
                   / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond6413]
            cond6414 = (
                (cond6412 | cond6413)
                & ((thX - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond6414] = (thX - self.soil.th_fc_adj[comp,:])[cond6414]

            # Update water content
            thnew[comp,:][cond641] = (thX - dthdt)[cond641]

            # Update cumulative drainage (mm), restrict to saturated hydraulic
            # conductivity and adjust excess drainage flow
            drainsum[cond641] += (dthdt * 1000 * self.soil.dz[comp])[cond641]
            cond6415 = (cond641 & (drainsum > ksat[comp,:]))
            excess[cond6415] += (drainsum - ksat[comp,:])[cond6415]
            drainsum[cond6415] = ksat[comp,:][cond6415]

            # ##################################################################

            cond642 = (
                cond64
                & np.logical_not(cond641)
                & (thnew[comp,:] > self.soil.th_fc_adj[comp,:]))
            
            # Calculate drainage ability for updated water content
            cond6421 = (cond642
                        & (thnew[comp,:] <= self.soil.th_fc_adj[comp,:]))
            dthdt[cond6421] = 0
            cond6422 = (
                cond642
                & np.logical_not(cond6421)
                & (thnew[comp,:] >= th_s[comp,:]))
            dthdt[cond6422] = (
                tau[comp,:]
                * (th_s[comp,:] - th_fc[comp,:]))[cond6422]
            cond6423 = (cond642 & np.logical_not(cond6421 | cond6422))
            dthdt[cond6423] = (
                tau[comp,:]
                * (th_s[comp,:] - th_fc[comp,:])
                * ((np.exp(thnew[comp,:] - th_fc[comp,:]) - 1)
                   / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond6423]
            cond6424 = (
                (cond6422 | cond6423)
                & ((thnew[comp,:] - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond6424] = (thnew[comp,:]
                               - self.soil.th_fc_adj[comp,:])[cond6424]
            
            # Update water content
            thnew[comp,:][cond642] = (thnew[comp,:] - dthdt)[cond642]

            # Update cumulative drainage (mm), restrict to saturated hydraulic
            # conductivity and adjust excess drainage flow
            drainsum[cond642] = (dthdt * 1000 * self.soil.dz[comp])[cond642]
            cond6425 = (cond642 & (drainsum > ksat[comp,:]))
            excess[cond6425] += (drainsum - ksat[comp,:])[cond6425]
            drainsum[cond6425] = ksat[comp,:][cond6425]

            # ##################################################################

            # Increase water content in compartment ii with cumulative drainage
            # from above
            cond65 = (np.logical_not(drainability) & (thX > th_s[comp,:]))
            thnew[comp,:][cond65] = (
                self.th[comp,:]
                + (drainsum / (1000 * self.soil.dz[comp])))[cond65]

            # Check new water content against hydraulic properties of soil layer
            cond651 = (cond65 & (thnew[comp,:] <= th_s[comp,:]))

            # "Calculate new drainage ability
            cond6511 = (cond651 & (thnew[comp,:] > self.soil.th_fc_adj[comp,:]))

            cond65111 = (cond6511 & (thnew[comp,:] <= th_fc[comp,:]))
            dthdt[cond65111] = 0

            cond65112 = (
                cond6511
                & np.logical_not(cond65111)
                & (thnew[comp,:] >= th_s[comp,:]))
            dthdt[cond65112] = (
                tau[comp,:]
                * (th_s[comp,:] - th_fc[comp,:]))[cond65112]

            cond65113 = (cond6511 & np.logical_not(cond65111 | cond65112))
            dthdt[cond65113] = (
                tau[comp,:]
                * (th_s[comp,:] - th_fc[comp,:])
                * ((np.exp(thnew[comp,:] - th_fc[comp,:]) - 1)
                   / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond65113]

            cond65114 = (
                (cond65112 | cond65113)
                & ((thnew[comp,:] - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond65114] = (thnew[comp,:] -
                                self.soil.th_fc_adj[comp,:])[cond65114]

            # Update water content
            thnew[comp,:][cond6511] -= (dthdt)[cond6511]

            # Update cumulative drainage (mm), restrict to saturated hydraulic
            # conductivity and adjust excess drainage flow
            drainsum[cond6511] = (dthdt * 1000 * self.soil.dz[comp])[cond6511]
            cond65115 = (cond6511 & (drainsum > ksat[comp,:]))
            excess[cond65115] += (drainsum - ksat[comp,:])[cond65115]
            drainsum[cond65115] = ksat[comp,:][cond65115]

            # ##################################################################

            cond652 = (cond65 & (thnew[comp,:] > th_s[comp,:]))

            # Calculate excess drainage above saturation
            excess[cond652] = (
                (thnew[comp,:] - th_s[comp,:])
                * 1000 * self.soil.dz[comp])[cond652]

            # "Calculate drainage ability for updated water content"
            cond6521 = (cond652
                        & (thnew[comp,:] <= self.soil.th_fc_adj[comp,:]))
            dthdt[cond6521] = 0

            cond6522 = (
                cond652
                & (np.logical_not(cond6521) & (thnew[comp,:] >= th_s[comp,:])))
            dthdt[cond6522] = (
                tau[comp,:]
                * (th_s[comp,:] - th_fc[comp,:]))[cond6522]
            cond6523 = (cond652 & np.logical_not(cond6521 | cond6522))
            dthdt[cond6523] = (
                tau[comp,:]
                * (th_s[comp,:] - th_fc[comp,:])
                * ((np.exp(thnew[comp,:] - th_fc[comp,:]) - 1)
                   / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond6523]

            cond6524 = (
                (cond6522 | cond6523)
                & ((thnew[comp,:] - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond6524] = (thnew[comp,:]
                             - self.soil.th_fc_adj[comp,:])[cond6524]

            # Update water content
            thnew[comp,:][cond652] = (th_s[comp,:] - dthdt)[cond652]

            # Update drainage, maximum drainage, excess drainage
            draincomp[cond652] = (dthdt * 1000 * self.soil.dz[comp])[cond652]
            drainmax[cond652] = (dthdt * 1000 * prethick)[cond652]
            drainmax[cond652] = np.clip(drainmax, None, excess)[cond652]
            excess[cond652] -= drainmax[cond652]

            # Update cumulative drainage (mm), restrict to saturated hydraulic
            # conductivity and adjust excess drainage flow
            drainsum[cond652] = (draincomp + drainmax)[cond652]
            cond6525 = (cond652 & (drainsum > ksat[comp,:]))
            excess[cond6525] += (drainsum - ksat[comp,:])[cond6525]
            drainsum[cond6525] = ksat[comp,:][cond6525]

            # ##################################################################

            # Store output flux from compartment ii
            self.FluxOut[comp,:] = drainsum

            # Redistribute excess in compartment above
            precomp = comp + 1
            while (np.any(excess > 0)) & (precomp != 0):

                # Include condition here so that it is updated
                cond7 = (excess > 0)
                
                # Update compartment counter
                precomp -= 1

                # Update flux from compartment
                if (precomp < comp):
                    self.FluxOut[precomp,:][cond7] -= excess[cond7]

                # Increase water content to store excess
                thnew[precomp,:][cond7] += (
                    excess / (1000 * self.soil.dz[precomp]))[cond7]

                # Limit water content to saturation and adjust excess counter
                cond71 = (cond7 & (thnew[precomp,:] > th_s[precomp,:]))
                excess[cond71] = (
                    (thnew[precomp,:] - th_s[precomp,:])
                    * 1000 * self.soil.dz[precomp])[cond71]                
                thnew[precomp,:][cond71] = th_s[precomp,:][cond71]
                
                cond72 = (cond7 & np.logical_not(cond71))
                excess[cond72] = 0


        # Update state variables
        self.DeepPerc = drainsum
        self.th = thnew

    def rainfall_partition(self, rotation, meteo):
        """Function to partition rainfall into surface runoff and 
        infiltration using the curve number approach.
        """
        # Add rotation dimension to precipitation
        P = meteo.precipitation[None,:,:] * np.ones((self.nRotation))[:,None,None]
        
        # Expand soil properties to compartments
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]

        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = self.soil.dz[:,None,None,None] * arr_ones
        dzsum = self.soil.dzsum[:,None,None,None] * arr_ones
        zcn = self.soil.zCN[None,:,:,:] * np.ones((self.nComp))[:,None,None,None]

        # Calculate runoff
        cond1 = ((rotation.Bunds == 0) | (rotation.zBund < 0.001))
        self.DaySubmerged[cond1] = 0
        cond11 = (cond1 & (self.soil.AdjCN == 1))

        # Check which compartment cover depth of top soil used to adjust
        # curve number
        comp_sto = ((dzsum - dz) < zcn)
        cond111 = np.all((comp_sto == False), axis=0)
        cond111 = np.broadcast_to(cond111, comp_sto.shape)
        comp_sto[cond111] = True

        # Calulcate weighting factors by compartment
        dzsum[dzsum > zcn] = zcn[dzsum > zcn]
        wx = (1.016 * (1 - np.exp(-4.16 * (dzsum / zcn))))

        # xx is wx for the overlying layer, with the top layer equal to zero
        xx = np.concatenate((np.zeros((1, self.nRotation, self.nLat, self.nLon)), wx[:-1,:]), axis=0)
        wrel = np.clip((wx - xx), 0, 1)

        # Calculate relative wetness of topsoil
        th = np.maximum(th_wp, self.th)

        # Multiply by comp_sto to ensure that compartments not used for
        # curve number adjustment are set to zero
        wet_top_comp = wrel * ((th - th_wp) / (th_fc - th_wp)) * comp_sto
        wet_top = np.sum(wet_top_comp, axis=0)
        wet_top = np.clip(wet_top, 0, 1)

        # Make adjustment to curve number
        CN = self.soil.CN
        CN[cond11] = np.round(self.soil.CNbot
                              + (self.soil.CNtop - self.soil.CNbot)
                              * wet_top)[cond11]

        # Partition rainfall into runoff and infiltration
        S = (25400 / CN) - 254
        term = P - ((5 / 100) * S)
        
        cond12 = (cond1 & (term <= 0))
        self.Runoff[cond12] = 0
        self.Infl[cond12] = P[cond12]

        cond13 = (cond1 & np.logical_not(cond12))
        self.Runoff[cond13] = ((term ** 2) / (P + (1 - (5 / 100)) * S))[cond13]
        self.Infl[cond13] = (P - self.Runoff)[cond13]

        # If there are bunds on the field then there is no runoff
        cond2 = np.logical_not(cond1)
        self.Runoff[cond2] = 0
        self.Infl[cond2] = P[cond2]

    def root_zone_water(self, soil_water):
        """Function to calculate actual and total available water in the 
        root zone at current time step
        """
        # Expand soil properties to compartments
        th_s   = self.soil.th_s[self.soil.layerIndex,:]
        th_fc  = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp  = self.soil.th_wp[self.soil.layerIndex,:]
        th_dry = self.soil.th_dry[self.soil.layerIndex,:]

        # Add rotation, lat, lon dimensions to dz and dzsum
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = self.soil.dz[:,None,None,None] * arr_ones
        dzsum = self.soil.dzsum[:,None,None,None] * arr_ones

        # Calculate root zone water content and available water
        rootdepth = np.maximum(self.Zmin, self.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        comp_sto = ((dzsum - dz) < rootdepth)
        
        # Fraction of compartment covered by root zone (zero in compartments
        # NOT covered by the root zone)
        factor = np.clip(rootdepth / dzsum, 0, 1)
        factor[np.logical_not(comp_sto)] = 0

        # Water storages in root zone (mm) - initially compute value in each
        # compartment, then sum to get overall root zone storages
        Wr_comp = factor * 1000 * soil_water.th * dz
        WrS_comp = factor * 1000 * th_s * dz
        WrFC_comp = factor * 1000 * th_fc * dz
        WrWP_comp = factor * 1000 * th_wp * dz
        WrDry_comp = factor * 1000 * th_dry * dz

        # Water storage in root zone at aeration stress threshold (mm)
        WrAer_comp = factor * 1000 * (th_s - (self.Aer / 100)) * dz

        Wr = np.sum(Wr_comp, axis=0)
        Wr[Wr < 0] = 0
        WrS = np.sum(WrS_comp, axis=0)
        WrFC = np.sum(WrFC_comp, axis=0)
        WrWP = np.sum(WrWP_comp, axis=0)
        WrDry = np.sum(WrDry_comp, axis=0)
        WrAer = np.sum(WrAer_comp, axis=0)

        # Convert depths to m3/m3
        self.thRZ['Act'] = Wr / (rootdepth * 1000)
        self.thRZ['Sat'] = WrS / (rootdepth * 1000)
        self.thRZ['Fc'] = WrFC / (rootdepth * 1000)
        self.thRZ['Wp'] = WrWP / (rootdepth * 1000)
        self.thRZ['Dry'] = WrDry / (rootdepth * 1000)
        self.thRZ['Aer'] = WrAer / (rootdepth * 1000)

        # Calculate total available water and root zone depletion
        self.TAW = np.clip((WrFC - WrWP), 0, None)
        self.Dr = np.clip((WrFC - Wr), 0, None)

    def irrigation(self, rotation, meteo):
        """Function to get irrigation depth for the current day"""

        # Expand soil properties to compartments
        th_s   = self.soil.th_s[self.soil.layerIndex,:]
        th_fc  = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp  = self.soil.th_wp[self.soil.layerIndex,:]
        th_dry = self.soil.th_dry[self.soil.layerIndex,:]

        SMT = np.concatenate((self.SMT1[None,:],
                              self.SMT2[None,:],
                              self.SMT3[None,:],
                              self.SMT4[None,:]), axis=0)

        # Calculate root zone water content and depletion
        self.root_zone_water()

        # Determine adjustment for inflows and outflows on current day
        cond1 = (self.thRZ['Act'] > self.thRZ['Fc'])
        rootdepth = np.zeros((self.nRotation, self.nLat, self.nLon))
        rootdepth[cond1] = np.maximum(self.Zmin, self.Zroot)[cond1]
        AbvFc = np.zeros((self.nRotation, self.nLat, self.nLon))
        AbvFc[cond1] = ((self.thRZ['Act'] - self.thRZ['Fc']) * 1000 * rootdepth)[cond1]
        WCadj = self.Tpot + self.Epot - meteo.precipitation + self.Runoff - AbvFc

        # Determine irrigation depth (mm/day) to be applied

        # Update growth stage if it is first day of a growing season
        cond2 = (growing_season_index & (self.DAP == 1))
        self.GrowthStage[cond2] = 1

        # Run irrigation depth calculation
        condX = (growing_season_index & (self.IrrMethod == 0))
        self.Irr[condX] = 0
        
        # If irrigation is based on soil moisture, get the soil moisture
        # target for the current growth stage and determine threshold to
        # initiate irrigation
        cond3 = (growing_season_index & (self.IrrMethod == 1))
        I,J,K = np.ogrid[:self.nRotation,:self.nLat,:self.nLon]
        growth_stage_index = self.GrowthStage - 1
        growth_stage_index = growth_stage_index.astype(int)
        SMT = SMT[growth_stage_index,I,J,K]
        IrrThr = (1 - SMT / 100) * self.TAW
                
        # If irrigation is based on a fixed interval, get number of days in
        # growing season so far (subtract 1 so that we always irrigate first
        # on day 1 of each growing season)
        cond4 = (growing_season_index & (self.IrrMethod == 2))
        nDays = self.DAP - 1

        # combine irrigation methods 1 and 2
        cond5 = (cond3 | cond4)
        
        # Working on a copy, adjust depletion for inflows and outflows - same
        # for both soil moisture and interval based irrigation
        Dr = self.Dr
        Dr[cond5] += WCadj[cond5]
        Dr[cond5] = np.clip(Dr, 0, None)[cond5]

        # check if conditions for irrigation method 1 or 2 are met
        cond6 = ((cond3 & (Dr > IrrThr))
                 | (cond4 & ((nDays % IrrInterval) == 0)))
        
        IrrReq = Dr
        EffAdj = ((100 - self.AppEff) + 100) / 100
        IrrReq *= EffAdj
        self.Irr[cond6] = np.clip(IrrReq, 0, self.MaxIrr)[cond6]

        # If irrigation is based on a pre-defined schedule then the irrigation
        # requirement for each rotation is read from a netCDF file. Note that if
        # the option 'irrScheduleFileNC' is None, then nothing will be imported
        # and the irrigation requirement will be set to zero
        cond8 = (growing_season_index & (self.IrrMethod == 3))
        if self.irrScheduleFileNC != None:
            IrrReq = vos.netcdf2PCRobjClone(self.irrScheduleFileNC,\
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
        cond9 = (growing_season_index & (self.IrrMethod == 4))
        self.Irr[cond9] = 0

        # No irrigation outside the growing season
        self.Irr[np.logical_not(growing_season_index)] = 0

        # Update cumulative irrigation counter for growing season
        self.IrrCum += self.Irr
        
    def infiltration(self, rotation):
        """Function to infiltrate incoming water (rainfall and 
        irrigation)
        """
        # Expand soil properties to compartments
        ksat   = self.soil.ksat[self.soil.layerIndex,:]
        th_s   = self.soil.th_s[self.soil.layerIndex,:]
        th_fc  = self.soil.th_fc[self.soil.layerIndex,:]
        tau    = self.soil.tau[self.soil.layerIndex,:]

        # initialise variables
        ToStore = np.zeros((self.nRotation, self.nLat, self.nLon))
        RunoffIni = np.zeros((self.nRotation, self.nLat, self.nLon))

        # Work on a copy of current water content
        thnew = self.th
        
        # Update infiltration rate for irrigation
        # Note: irrigation amount adjusted for specified application efficiency
        self.Infl += rotation.Irr * (rotation.AppEff / 100)  # TODO: has efficiency not already been accounted for?

        # Determine surface storage if bunds are present
        cond1 = (rotation.Bunds == 1)
        cond11 = (cond1 & (rotation.zBund > 0.001))
        InflTot = self.Infl + self.SurfaceStorage

        # Update surface storage and infiltration storage
        cond111 = (cond11 & (InflTot > 0))

        # Infiltration limited by saturated hydraulic conductivity of surface
        # soil layer; additional water ponds on surface
        cond1111 = (cond111 & (InflTot > ksat[0,:]))
        ToStore[cond1111] = ksat[0,:][cond1111]
        self.SurfaceStorage[cond1111] = (InflTot - ksat[0,:])[cond1111]

        # Otherwise all water infiltrates and surface storage becomes zero
        cond1112 = (cond111 & np.logical_not(cond1111))
        ToStore[cond1112] = InflTot[cond1112]
        self.SurfaceStorage[cond1112] = 0
        
        # Calculate additional runoff if water overtops bunds
        cond1113 = (cond111 & (self.SurfaceStorage > (rotation.zBund * 1000)))
        RunoffIni[cond1113] = (self.SurfaceStorage
                               - (rotation.zBund * 1000))[cond1113]
        self.SurfaceStorage[cond1113] = (rotation.zBund * 1000)[cond1113]

        # Otherwise excess water does not overtop bunds and there is no runoff
        cond1114 = (cond111 & np.logical_not(cond1113))
        RunoffIni[cond1114] = 0
        
        # If bunds are not on field then infiltration is limited by saturated
        # hydraulic conductivity of top soil layer
        cond2 = (rotation.Bunds == 0)
        cond21 = (cond2 & (self.Infl > ksat[0,:]))
        ToStore[cond21] = ksat[0,:][cond21]
        RunoffIni[cond21] = (self.Infl - ksat[0,:])[cond21]

        # Otherwise all water infiltrates
        cond22 = (cond2 & np.logical_not(cond21))
        ToStore[cond22] = self.Infl[cond22]
        RunoffIni[cond22] = 0

        # Initialize counters
        comp = 0
        DeepPerc = np.zeros((self.nRotation, self.nLat, self.nLon))
        Runoff = np.zeros((self.nRotation, self.nLat, self.nLon))

        cond3_ini = (ToStore > 0)

        while (np.any(ToStore > 0) & (comp < self.soil.nComp)):

            # Update condition
            cond3 = (ToStore > 0)
            
            # Calculate saturated drainage ability, drainage factor and
            # required drainage ability
            dthdtS = tau[comp,:] * (th_s[comp,:] - th_fc[comp,:])
            factor = ksat[comp,:] / (dthdtS * 1000 * self.soil.dz[comp])
            dthdt0 = ToStore / (1000 * self.soil.dz[comp])

            # Initialize water content for current layer
            theta0 = np.zeros((self.nRotation, self.nLat, self.nLon))  # initialise

            # "Check drainage ability
            cond31 = (cond3 & (dthdt0 < dthdtS))

            # Calculate water content, thX, needed to meet drainage dthdt0
            cond311 = (cond31 & (dthdt0 < dthdtS))
            theta0[cond311] = self.th_fc_adj[comp,:][cond311]
            cond312 = (cond31 & np.logical_not(cond311))
            A = (1 + ((dthdt0 * (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1))
                      / (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]))))
            theta0[cond312] = (self.soil.th_fc[comp,:] + np.log(A))[cond312]

            # Limit thX to between saturation and field capacity
            cond313 = (cond31 & (theta0 > th_s[comp,:]))
            theta0[cond313] = th_s[comp,:][cond313]
            cond314 = (cond31 & (theta0 < self.soil.th_fc_adj[comp,:]))
            theta0[cond314] = self.soil.th_fc_adj[comp,:][cond314]
            dthdt0[cond314] = 0

            # Limit water content and drainage to saturation
            cond32 = (cond3 & np.logical_not(cond31))
            theta0[cond32] = th_s[comp,:][cond32]
            dthdt0[cond32] = dthdtS[cond32]

            # Calculate maximum water flow through compartment and total
            # drainage
            drainmax = factor * dthdt0 * 1000 * self.soil.dz[comp]
            drainage = drainmax + self.FluxOut[comp,:]

            # Limit drainage to saturated hydraulic conductivity
            cond33 = (cond3 & (drainage > ksat[comp,:]))
            drainmax[cond33] = (ksat[comp,:] - self.FluxOut[comp,:])[cond33]

            # Calculate difference between threshold and current water contents
            diff = theta0 - self.th[comp,:]

            cond34 = (cond3 & (diff > 0))
            thnew[comp,:][cond34] += (ToStore / (1000 * self.soil.dz[comp]))[cond34]

            # Remaining water that can infiltrate to compartments below
            cond341 = (cond34 & (thnew[comp,:] > theta0))
            ToStore[cond341] = (
                (thnew[comp,:] - theta0)
                * 1000
                * self.soil.dz[comp])[cond341]
            thnew[comp,:][cond341] = theta0[cond341]

            # Otherwise all infiltrating water has been stored
            cond342 = (cond34 & np.logical_not(cond341))
            ToStore[cond342] = 0

            # Update outflow from current compartment (drainage + infiltration
            # flows)
            self.FluxOut[comp,:][cond3] += ToStore[cond3]

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
                self.FluxOut[precomp,:][cond35] = (
                    self.FluxOut[precomp,:] - excess)[cond35]

                # Update water content and limit to saturation
                thnew[precomp,:][cond35] += (
                    excess / (1000 * self.soil.dz[precomp]))[cond35]
                cond351 = (cond35 & (thnew[precomp,:] > th_s[precomp,:]))
                excess[cond351] = ((thnew[precomp,:] - th_s[precomp,:]) * 1000 * self.soil.dz[precomp])[cond351]
                thnew[precomp,:][cond351] = th_s[precomp,:][cond351]
                cond352 = (cond35 & np.logical_not(cond351))
                excess[cond352] = 0

            # Any leftover water not stored becomes runoff
            cond36 = (cond3 & (excess > 0))
            Runoff[cond36] += excess[cond36]
        
        # Infiltration left to store after bottom compartment becomes deep
        # percolation (mm)
        DeepPerc = ToStore

        # Otherwise if ToStore equals zero there is no infiltration
        cond4 = np.logical_not(cond3_ini)
        Runoff[cond4] = 0
        DeepPerc[cond4] = 0
        
        # Update total runoff
        Runoff += RunoffIni

        # Update surface storage (if bunds are present)
        cond5 = ((Runoff > RunoffIni)
                 & (self.Bunds == 1)
                 & rotation.zBund > 0.001)
        self.SurfaceStorage += (Runoff - RunoffIni)

        # Limit surface storage to bund height: additional water above top of
        # bunds becomes runoff, and surface storage equals bund height
        cond51 = (cond5 & (self.SurfaceStorage > (rotation.zBund * 1000)))
        Runoff[cond51] = (
            RunoffIni
            + (self.SurfaceStorage - (rotation.zBund * 1000)))[cond51]
        self.SurfaceStorage[cond51] = (rotation.zBund * 1000)[cond51]
        cond52 = (cond5 & np.logical_not(cond51))
        Runoff[cond52] = RunoffIni[cond52]

        # Store updated water contents
        self.th = thnew

        # Update deep percolation, surface runoff, and infiltration values
        self.DeepPerc += DeepPerc
        self.Infl -= Runoff
        self.Runoff += Runoff

    def capillary_rise(self, groundwater):
        """Function to calculate capillary rise from a shallow 
        groundwater table
        """
        # Expand soil properties to compartments
        ksat  = self.soil.ksat[self.soil.layerIndex,:]
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]
        aCR   = self.soil.aCR[self.soil.layerIndex,:]
        bCR   = self.soil.bCR[self.soil.layerIndex,:]

        # Add rotation, lat, lon dimensions to dz and dzsum
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = self.soil.dz[:,None,None,None] * arr_ones
        dzsum = self.soil.dzsum[:,None,None,None] * arr_ones
        
        if groundwater.WaterTable:

            # Get maximum capillary rise for bottom compartment
            zBot = np.sum(self.soil.dz)
            zBotMid = zBot - (self.soil.dz[-1] / 2)

            MaxCR = np.zeros((self.nRotation, self.nLat, self.nLon))
            cond1 = ((ksat[-1,:] > 0) & (self.zGW > 0) & ((self.zGW - zBotMid) < 4))
            cond11 = (cond1 & (zBotMid >= self.zGW))
            MaxCR[cond11] = 99            
            cond12 = (cond1 & np.logical_not(cond11))
            MaxCR[cond12] = (np.exp((np.log(self.zGW - zBotMid) - bCR[-1,:]) / aCR[-1,:]))[cond12]
            MaxCR = np.clip(MaxCR, None, 99)

            # Find top of next soil layer that is not within modelled soil
            # profile
            # Note: self.soil.layerIndex contains the index of the layer
            # corresponding to each compartment; hence the last value is the
            # layer corresponding to the bottom compartment.
            zTopLayer = np.cumsum(self.soil.zLayer[np.arange(0, (self.soil.layerIndex[-1] + 1))])

            # Check for restrictions on upward flow caused by properties of
            # compartments that are not modelled in the soil water balance

            # Layer number corresponding to bottom compartment
            layeri = self.soil.layerIndex[-1]
            LimCR = np.zeros((nr, nlat, nlon))

            while np.any(zTopLayer < self.zGW) & (layeri < (self.soil.nLayer - 1)):
                layeri += 1     # layer fully below bottom compartment
                cond2 = (zTopLayer < self.zGW)
                cond21 = (
                    cond2
                    & (self.soil.ksat[layeri,:] > 0)
                    & (self.zGW > 0)
                    & ((self.zGW - zTopLayer) < 4))
                cond211 = (cond21 & (zTopLayer >= self.zGW))
                LimCR[cond211] = 99
                cond212 = (cond21 & np.logical_not(cond211))
                LimCR[cond212] = (
                    np.exp((np.log(self.zGW - zTopLayer)
                            - self.soil.bCR[layeri,:])
                           / self.soil.aCR[layeri,:]))[cond7]
                LimCR = np.clip(LimCR, None, 99)
                MaxCR = np.clip(MaxCR, None, LimCR)
                zTopLayer += self.soil.zLayer[layeri]

            # Calculate capillary rise, starting at the bottom compartment
            compi = self.soil.nComp - 1

            # Initialize capillary rise counter
            WCr = np.zeros((self.nRotation, self.nLat, self.nLon))
            while ((np.any(np.round(MaxCR * 1000) > 0))
                   & (np.any(np.round(self.FluxOut[compi,:] * 1000) == 0))
                   & (compi >= 0)):
                
                # Proceed upwards until maximum capillary rise occurs, soil
                # surface is reached, or encounter a compartment where downward
                # drainage/infiltration has already occurred on current day
                cond3 = ((np.round(MaxCR * 1000) > 0)
                         & (np.round(self.FluxOut[compi,:] * 1000) == 0))
                Df = np.zeros((nr, nlat, nlon))
                cond31 = (cond3 & ((self.th[compi,:] >= th_wp[compi,:]) & (self.soil.fshape_cr > 0)))
                Df[cond31] = (
                    1 - (((self.th[compi,:] - th_wp[compi,:])
                          / (self.soil.th_fc_adj[compi,:] - th_wp[compi,:]))
                         ** self.soil.fshape_cr))[cond31]
                Df = np.clip(Df, 0, 1)
                cond32 = (cond3 & np.logical_not(cond31))
                Df[cond32] = 1
                
                # Calculate relative hydraulic conductivity
                thThr = (th_wp[compi,:] + th_fc[compi,:]) / 2
                Krel = np.zeros((self.nRotation, self.nLat, self.nLon))
                cond33 = (cond3 & (self.th[compi,:] < thThr))                
                cond331 = (cond33 & ((self.th[compi,:] <= th_wp[compi,:]) | (thThr <= th_wp[compi,:])))
                Krel[cond331] = 0
                cond332 = (cond33 & np.logical_not(cond331))
                Krel[cond332] = (
                    (self.th[compi,:] - th_wp[compi,:])
                    / (thThr - th_wp[compi,:]))[cond332]
                
                cond34 = (cond3 & np.logical_not(cond33))
                Krel[cond34] = 1

                # Check if room is available to store water from capillary rise
                dth = np.zeros((self.nRotation, self.nLat, self.nLon))
                dth[cond3] = (
                    self.soil.th_fc_adj[compi,:]
                    - self.th[compi,:])[cond3]
                
                dth = np.round((dth * 1000) / 1000)

                # Store water if room is available
                dthMax = np.zeros((self.nRotation, self.nLat, self.nLon))
                CRcomp = np.zeros((self.nRotation, self.nLat, self.nLon))
                cond35 = (
                    cond3
                    & ((dth > 0) & ((zBot - (self.soil.dz[-1] / 2)) < self.zGW)))
                
                dthMax[cond35] = (Krel * Df * MaxCR
                                  / (1000 * self.soil.dz[compi]))[cond35]
                cond351 = (cond35 & (dth >= dthMax))
                self.th[compi,:][cond351] += dthMax[cond351]
                CRcomp[cond351] = (dthMax * 1000 * self.soil.dz[compi])[cond351]
                MaxCR[cond351] = 0
                cond352 = (cond35 & np.logical_not(cond351))
                self.th[compi,:][cond352] = self.soil.th_fc_adj[compi,:][cond352]
                CRcomp[cond352] = (dth * 1000 * self.soil.dz[compi])[cond352]
                MaxCR[cond352] = ((Krel * MaxCR) - CRcomp)[cond352]
                WCr[cond35] += CRcomp[cond35]

                # Update bottom elevation of compartment and compartment counter
                zBot -= self.soil.dz[compi]
                compi -= 1

                # Update restriction on maximum capillary rise
                if compi >= 0:
                    zBotMid = zBot - (self.soil.dz[compi] / 2)
                    cond36 = (
                        cond3 & ((ksat[compi,:] > 0)
                                 & (self.zGW > 0)
                                 & ((self.zGW - zBotMid) < 4)))
                    
                    cond361 = (cond36 & (zBotMid >= self.zGW))
                    LimCR[cond361] = 99
                    cond362 = (cond36 & np.logical_not(cond361))
                    LimCR[cond362] = (
                        np.exp((np.log(self.zGW - zBotMid) - bCR[compi,:])
                               / aCR[compi,:]))[cond362]

                    LimCR = np.clip(LimCR, None, 99)
                    cond37 = (cond3 & np.logical_not(cond36))
                    LimCR[cond37] = 0
                    cond38 = (cond3 & (MaxCR > LimCR))
                    MaxCR[cond38] = LimCR[cond38]
            
            # Store total depth of capillary rise
            self.CrTot = WCr

    def germination(self):
        """Function to check if crop has germinated"""

        # Expand soil properties to compartments
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]

        # Add rotation, lat, lon dimensions to dz and dzsum
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))
        dz = self.soil.dz[:,None,None,None] * arr_ones
        dzsum = self.soil.dzsum[:,None,None,None] * arr_ones

        # Add compartment dimension to zGerm
        zgerm = self.soil.zGerm[None,:,:,:] * np.ones((self.nComp))[:,None,None,None]

        # Here we force zGerm to have a maximum value equal to the depth of the
        # deepest soil compartment
        zgerm[zgerm > np.sum(self.soil.dz)] = np.sum(self.soil.dz)
        
        # Find compartments covered by top soil layer affecting germination
        comp_sto = ((dzsum - dz) < zgerm)

        # Calculate water content in top soil layer
        arr_zeros = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))
        Wr_comp = arr_zeros
        WrFC_comp = arr_zeros
        WrWP_comp = arr_zeros

        # Determine fraction of compartment covered by top soil layer
        factor = 1 - ((dzsum - zgerm) / dz) 
        factor = np.clip(factor, 0, 1) * growing_season_index * comp_sto

        # Increment water storages (mm)
        Wr_comp = (factor * 1000 * self.th * dz)
        Wr_comp = np.clip(Wr_comp, 0, None)
        Wr = np.sum(Wr_comp, axis=0)
        WrFC_comp = (factor * 1000 * th_fc * dz)
        WrFC = np.sum(WrFC_comp, axis=0)
        WrWP_comp = (factor * 1000 * th_wp * dz)
        WrWP = np.sum(WrWP_comp, axis=0)

        # Calculate proportional water content
        WcProp = 1 - ((WrFC - Wr) / (WrFC - WrWP))

        # Check if water content is above germination threshold
        cond4 = (growing_season_index
                 & (WcProp >= self.GermThr)
                 & (np.logical_not(self.Germination)))
        self.Germination[cond4] = True

        # Increment delayed growth time counters if germination is yet to occur
        cond5 = (growing_season_index & (np.logical_not(self.Germination)))
        self.DelayedCDs[cond5] += 1
        self.DelayedGDDs[cond5] += self.GDD[cond5]
                 
        # Not in growing season so no germination calculation is performed
        self.Germination[np.logical_not(growing_season_index)] = False
        self.DelayedCDs[np.logical_not(growing_season_index)] = 0
        self.DelayedGDDs[np.logical_not(growing_season_index)] = 0

    def growth_stage(self):
        """Function to calculate number of growing degree days on 
        current day
        """
        if self.crop.CalendarType == 1:
            tAdj = self.DAP - self.DelayedCDs
        elif self.crop.CalendarType == 2:
            tAdj = self.GDDcum - self.DelayedGDDs

        # Update growth stage
        cond1 = (growing_season_index & (tAdj <= self.Canopy10Pct))
        cond2 = (growing_season_index & (tAdj <= self.MaxCanopy))
        cond3 = (growing_season_index & (tAdj <= self.Senescence))
        cond4 = (growing_season_index & (tAdj > self.Senescence))

        self.GrowthStage[cond1] = 1
        self.GrowthStage[cond2] = 2
        self.GrowthStage[cond3] = 3
        self.GrowthStage[cond4] = 4

    def root_development(self, groundwater):
        """Function to calculate root zone expansion"""

        cond1 = growing_season_index
        
        # If today is the first day of season, root depth is equal to minimum depth
        cond1 = (growing_season_index & (self.DAP == 1))
        self.Zroot[cond1] = self.Zmin[cond1]

        # Adjust time for any delayed development
        if self.CalendarType == 1:
            tAdj = (self.DAP - self.DelayedCDs)
        elif self.CalendarType == 2:
            tAdj = (self.GDDcum - self.DelayedGDDs)
            
        # Calculate root expansion
        Zini = self.Zmin * (self.PctZmin / 100)
        t0 = np.round(self.Emergence / 2)
        tmax = self.MaxRooting
        if self.CalendarType == 1:
            tOld = (tAdj - 1)
        elif self.CalendarType == 2:
            tOld = (tAdj - self.GDD)

        tAdj[np.logical_not(growing_season_index)] = 0
        tOld[np.logical_not(growing_season_index)] = 0

        # Potential root depth on previous day
        ZrOld = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond2 = (growing_season_index & (tOld >= tmax))
        ZrOld[cond2] = self.Zmax[cond2]
        cond3 = (growing_season_index & (tOld <= t0))
        ZrOld[cond3] = Zini[cond3]
        cond4 = (growing_season_index & (np.logical_not(cond2 | cond3)))
        X = (tOld - t0) / (tmax - t0)
        ZrOld[cond4] = ((Zini + (self.Zmax - Zini))
                        * X ** (1 / self.fshape_r))[cond4]
        
        cond5 = (growing_season_index & (ZrOld < self.Zmin))
        ZrOld[cond5] = self.Zmin[cond5]

        # Potential root depth on current day
        Zr = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond6 = (growing_season_index & (tAdj >= tmax))
        Zr[cond6] = self.Zmax[cond6]
        cond7 = (growing_season_index & (tAdj <= t0))
        Zr[cond7] = Zini[cond7]
        cond8 = (growing_season_index & (np.logical_not(cond6 | cond7)))
        ZrOld[cond8] = ((Zini + (self.Zmax - Zini))
                        * X ** (1 / self.fshape_r))[cond8]
        
        cond9 = (growing_season_index & (Zr < self.Zmin))
        Zr[cond9] = self.Zmin[cond9]
        
        # Determine rate of change, adjust for any stomatal water stress
        dZr = Zr - ZrOld
        cond10 = (growing_season_index & (self.TrRatio < 0.9999))
        cond101 = (cond10 & (self.fshape_ex >= 0))        
        dZr[cond101] = (dZr * self.TrRatio)[cond101]
        cond102 = (cond10 & np.logical_not(cond101))
        fAdj = ((np.exp(self.TrRatio * self.fshape_ex) - 1)
                / (np.exp(self.fshape_ex) - 1))
        dZr[cond102] = (dZr * fAdj)[cond102]

        # Adjust root expansion for failure to germinate (roots cannot expand
        # if crop has not germinated)
        dZr[np.logical_not(self.Germination)] = 0

        # Get new rooting depth
        self.Zroot[growing_season_index] = (self.Zroot + dZr)[growing_season_index]

        cond11 = (growing_season_index & (self.soil.zRes > 0))
        cond111 = (cond11 & (self.Zroot > self.soil.zRes))
        self.rCor[cond111] = (
            (2 * (self.Zroot / self.soil.zRes)
             * ((self.SxTop + self.SxBot) / 2)
             - self.SxTop) / self.SxBot)[cond111]
        
        self.Zroot[cond111] = self.soil.zRes[cond111]

        # Limit rooting depth if groundwater table is present (roots cannot
        # develop below the water table)
        zGW = groundwater.zGW[None,:,:] * np.ones((nr))[:,None,None]
        cond12 = (growing_season_index & groundwater.WaterTable & (zGW > 0))
        cond121 = (cond12 & (self.Zroot > zGW))
        self.Zroot[cond121] = zGW[cond121]
        cond1211 = (cond121 & (self.Zroot < self.Zmin))
        self.Zroot[cond1211] = self.Zmin[cond1211]

        # No root system outside of growing season
        self.Zroot[np.logical_not(growing_season_index)] = 0

    def water_stress(self, meteo, beta):
        """Function to calculate water stress coefficients"""        
        p_up = np.concatenate(
            (self.p_up1[None,:],
             self.p_up2[None,:],
             self.p_up3[None,:],
             self.p_up4[None,:]), axis=0)

        p_lo = np.concatenate(
            (self.p_lo1[None,:],
             self.p_lo2[None,:],
             self.p_lo3[None,:],
             self.p_lo4[None,:]), axis=0)
        
        fshape_w = np.concatenate(
            (self.fshape_w1[None,:],
             self.fshape_w2[None,:],
             self.fshape_w3[None,:],
             self.fshape_w4[None,:]), axis=0)
        
        # Adjust stress thresholds for Et0 on current day (don't do this for
        # pollination water stress coefficient)

        et0 = (meteo.referencePotET[None,:,:] * np.ones((self.nRotation))[:,None,None])
        cond1 = (self.ETadj == 1)
        for stress in range(3):
            p_up[stress,:][cond1] = (
                p_up[stress,:]
                + (0.04 * (5 - et0))
                * (np.log10(10 - 9 * p_up[stress,:])))[cond1]
            
            p_lo[stress,:][cond1] = (
                p_lo[stress,:]
                + (0.04 * (5 - et0))
                * (np.log10(10 - 9 * p_lo[stress,:])))[cond1]

        # Adjust senescence threshold if early senescence triggered
        if beta:
            cond2 = (self.tEarlySen > 0)
            p_up[2,:][cond2] = (p_up[2,:] * (1 - (self.beta / 100)))[cond2]

        # Limit adjusted values
        p_up = np.clip(p_up, 0, 1)
        p_lo = np.clip(p_lo, 0, 1)

        # Calculate relative depletion
        Drel = np.zeros((4, self.nRotation, self.nLat, self.nLon))
        # No water stress
        cond1 = (self.Dr <= (p_up * self.TAW))
        Drel[cond1] = 0
        
        # Partial water stress
        cond2 = (self.Dr >  (p_up * self.TAW)) & (self.Dr < (p_lo * self.TAW))
        Drel[cond2] = (1 - ((p_lo - (self.Dr / self.TAW)) / (p_lo - p_up)))[cond2]
        
        # Full water stress
        cond3 = (self.Dr >= (p_lo * self.TAW))
        Drel[cond3] = 1         

        # Calculate root zone stress coefficients
        idx = np.arange(0,3)
        Ks = (1 - ((np.exp(Drel[idx,:] * fshape_w[idx,:]) - 1)
                   / (np.exp(fshape_w[idx,:]) - 1)))

        # Water stress coefficients (leaf expansion, stomatal closure,
        # senescence, pollination failure)
        self.Ksw['Exp'] = Ks[0,:]
        self.Ksw['Sto'] = Ks[1,:]
        self.Ksw['Sen'] = Ks[2,:]
        self.Ksw['Pol'] = 1 - Drel[3,:]

        # Mean water stress coefficient for stomatal closure
        self.Ksw['StoLin'] = 1 - Drel[1,:]

    def canopy_cover_development(self, CC0, CCx, CGC, CDC, dt, Mode):
        """Function to calculate canopy cover development by end of the 
        current simulation day
        """
        if Mode == 'Growth':
            CC = (CC0 * np.exp(CGC * dt))
            cond1 = (CC > (CCx / 2))
            CC[cond1] = (CCx - 0.25 * (CCx / CC0) * CCx * np.exp(-CGC * dt))[cond1]
            CC = np.clip(CC, None, CCx)
        elif Mode == 'Decline':
            CC = np.zeros((CC0.shape))
            cond2 = (CCx >= 0.001)
            CC[cond2] = (CCx * (1 - 0.05 * (np.exp(dt * (CDC / CCx)) - 1)))[cond2]

        CC = np.clip(CC, 0, 1)
        return CC

    def canopy_cover_required_time(self, CCprev, CC0, CCx, CGC, CDC, dt, tSum, Mode):
        """Function to find required time to reach CC at end of previous 
        day, given current CGC or CDC
        """
        if Mode == 'CGC':
            CGCx = np.zeros((self.nRotation, self.nLat, self.nLon))
            cond1 = (CCprev <= (CCx / 2))
            CGCx[cond1] = ((np.log(CCprev / CC0)) / (tSum - dt))[cond1]
            cond2 = np.logical_not(cond1)
            CGCx[cond2] = (
                (np.log((0.25 * CCx * CCx / CC0) / (CCx - CCprev)))
                / (tSum - dt))[cond2]
            tReq = (tSum - dt) * (CGCx / CGC)
        elif Mode == 'CDC':
            tReq = ((np.log(1 + (1 - CCprev / CCx) / 0.05)) / (CDC / CCx))

        return tReq
                    
    def adjust_CCx(self, CCprev, CC0, CCx, CGC, CDC, dt, tSum):
        """Function to adjust CCx value for changes in CGC due to water 
        stress during the growing season
        """
        # Get time required to reach CC on previous day, then calculate
        # adjusted CCx
        tCCtmp = self.canopy_cover_required_time(CCprev, CC0, CCx, CGC,
                                                 CDC, dt, tSum, 'CGC')
        cond1 = (tCCtmp > 0)
        tCCtmp[cond1] += ((self.CanopyDevEnd - tSum) + dt)[cond1]
        CCxAdj = self.canopy_cover_development(CC0, CCx, CGC,
                                               CDC, tCCtmp, 'Growth')
        CCxAdj[np.logical_not(cond1)] = 0
        return CCxAdj

    def update_CCx_and_CDC(self, CCprev, CDC, CCx, dt):
        """Function to update CCx and CDC parameter values for 
        rewatering in late season of an early declining canopy
        """
        CCxAdj = CCprev / (1 - 0.05 * (np.exp(dt * (CDC / CCx)) - 1))
        CDCadj = CDC * (CCxAdj / CCx)
        return CCxAdj,CDCadj

    def canopy_cover(self, meteo):
        """Function to simulate canopy growth/decline"""

        # Preallocate some variables
        CCxAdj = np.zeros((self.nRotation, self.nLat, self.nLon))
        CDCadj = np.zeros((self.nRotation, self.nLat, self.nLon))
        CCsen = np.zeros((self.nRotation, self.nLat, self.nLon))

        # Store initial condition
        self.CCprev = self.CC
        
        # Calculate root zone water content, and determine if water stress is
        # occurring
        self.root_zone_water()
        self.water_stress(meteo, beta = True)

        # Get canopy cover growth over time
        if self.crop.CalendarType == 1:
            tCC = self.DAP
            dtCC = np.ones((self.nRotation, self.nLat, self.nLon))
            tCCadj = self.DAP - self.DelayedCDs
        elif self.crop.CalendarType == 2:
            tCC = self.GDDcum
            dtCC = self.GDD
            tCCadj = self.GDDcum - self.DelayedGDDs

        # ######################################################################
        # Canopy development (potential)

        # No canopy development before emergence/germination or after maturity
        cond1 = (growing_season_index
                 & ((tCC < self.Emergence) | (np.round(self.Emergence) > self.Maturity)))
        self.CC_NS[cond1] = 0

        # Canopy growth can occur
        cond2 = (growing_season_index & (tCC < self.CanopyDevEnd))

        # Very small initial CC as it is first day or due to senescence. In this
        # case assume no leaf expansion stress
        cond21 = (cond2 & (self.CC_NS <= self.CC0))
        self.CC_NS[cond21] = (self.CC0 * np.exp(self.CGC * dtCC))[cond21]
        
        # Canopy growing
        cond22 = (cond2 & np.logical_not(cond21))
        tmp_tCC = tCC - self.Emergence
        tmp_CC = self.canopy_cover_development(self.CC0, self.CCx, self.CGC, self.CDC, tmp_tCC, 'Growth')
        self.CC[cond22] = tmp_CC[cond22]

        # Update maximum canopy cover size in growing season
        self.CCxAct_NS[cond2] = self.CC_NS[cond2]
        
        # No more canopy growth is possible or canopy in decline
        cond3 = (growing_season_index & (tCC > self.CanopyDevEnd))
        # Set CCx for calculation of withered canopy effects
        self.CCxW_NS[cond3] = self.CCxAct_NS[cond3]
        
        # Mid-season stage - no canopy growth, so do not update CC_NS
        cond31 = (cond3 & (tCC < self.Senescence))
        self.CCxAct_NS[cond31] = self.CC_NS[cond31]

        # Late-season stage - canopy decline
        cond32 = (cond3 & np.logical_not(cond31))
        tmp_tCC = tCC - self.Emergence
        tmp_CC = self.canopy_cover_development(self.CC0, self.CCx, self.CGC, self.CDC, tmp_tCC, 'Decline')
        self.CC[cond32] = tmp_CC[cond32]

        # ######################################################################
        # Canopy development (actual)

        # No canopy development before emergence/germination or after maturity
        cond4 = (growing_season_index
                 & ((tCCadj < self.Emergence) | (np.round(tCCadj) > self.Maturity)))
        self.CC[cond4] = 0

        # Canopy growth can occur
        cond5 = (growing_season_index
                 & (tCCadj < CanopyDevEnd)
                 & np.logical_not(cond4))
        cond51 = (cond5 & (self.CCprev <= self.CC0adj))

        # Very small initial CC as it is first day or due to senescence. In
        # this case, assume no leaf expansion stress
        self.CC[cond51] = (self.CC0adj * np.exp(self.CGC * dtCC))[cond51]

        # Canopy growing
        cond52 = (cond5 & np.logical_not(cond51))

        # Canopy approaching maximum size
        cond521 = (cond52 & (self.CCprev >= (0.9799 * self.CCx)))
        tmp_tCC = tCC - self.Emergence
        tmp_CC = self.canopy_cover_development(self.CC0, self.CCx, self.CGC, self.CDC, tmp_tCC, 'Growth')
        self.CC[cond521] = tmp_CC[cond521]
        self.CC0adj[cond521] = self.CC0[cond521]

        # Adjust canopy growth coefficient for leaf expansion water stress
        # effects
        cond522 = (cond52 & np.logical_not(cond521))
        CGCadj = self.CGC * self.Ksw['Exp']

        # Adjust CCx for change in CGC
        cond5221 = (cond522 & (CGCadj > 0))
        tmp_CCxAdj = self.adjust_CCx(self.CCprev, self.CC0adj, self.CCx,
                                        CGCadj, self.CDC, dtCC, tCCadj)
        CCxAdj[cond5221] = tmp_CCxAdj[cond5221]

        cond52211 = (cond5221 & (CCxAdj > 0))

        # Approaching maximum canopy size
        cond522111 = (cond52211 & (np.abs(self.CCprev - self.CCx) < 0.00001))
        tmp_tCC = tCC - self.Emergence
        tmp_CC = self.canopy_cover_development(self.CC0, self.CCx, self.CGC, self.CDC, tmp_tCC, 'Growth')
        self.CC[cond522111] = tmp_CC[cond522111]

        # Determine time required to reach CC on previous day, given CGCadj
        # value
        cond522112 = (cond52211 & np.logical_not(cond522111))
        tReq = self.canopy_cover_required_time(self.CCprev, self.CC0adj, CCxAdj,
                                       CGCadj, self.CDC, dtCC, tCCadj, 'CGC')
        tmp_tCC = tReq + dtCC

        # Determine new canopy size
        cond5221121 = (cond522112 & (tmp_tCC > 0))
        tmp_CC = self.canopy_cover_development(self.CC0adj, CCxAdj, CGCadj,
                                        self.CDC, tmp_tCC, 'Growth')
        self.CC[cond5221121] = tmp_CC[cond5221121]

        # No canopy growth
        cond5222 = (cond522 & np.logical_not(cond5221))

        # Update CC0 if current canopy cover if less than initial canopy cover
        # size at planting
        cond52221 = (cond5222 & (self.CC < self.CC0adj))
        self.CC0adj[cond52221] = self.CC[cond52221]

        # "Update actual maximum canopy cover size during growing season"
        cond53 = (cond5 & (self.CC > self.CCxAct))
        self.CCxAct[cond53] = self.CC[cond53]
        
        # No more canopy growth is possible or canopy is in decline
        cond6 = (growing_season_index & (tCCadj > CanopyDevEnd))

        # Mid-season stage - no canopy growth: update actual maximum canopy
        # cover size during growing season only (i.e. do not update CC)
        cond61 = (cond6 & (tCCadj < self.Senescence))
        cond611 = (cond61 & (self.CC > self.CCxAct))
        self.CCxAct[cond611] = self.CC[cond611]

        # Late season stage - canopy decline: update canopy decline coefficient
        # for difference between actual and potential CCx, and determine new
        # canopy size
        cond62 = (cond6 & np.logical_not(cond61))
        CDCadj[cond62] = (self.CDC * (self.CCxAct / self.CCx))[cond62]
        tmp_tCC = tCCadj - self.Senescence
        tmp_CC = self.canopy_cover_development(self.CC0adj, self.CCxAct, self.CGC,
                                        CDCadj, tmp_tCC, 'Decline')
        self.CC[cond62] = tmp_CC[cond62]

        # Check for crop growth termination: if the following conditions are
        # met, the crop has died
        cond63 = (cond6 & ((self.CC < 0.001) & np.logical_not(self.CropDead)))
        self.CC[cond63] = 0
        self.CropDead[cond63] = True

        # ######################################################################
        # Canopy senescence due to water stress (actual)

        cond7 = (growing_season_index & (tCCadj >= self.Emergence))

        # Check for early canopy senescence starting/continuing due to severe
        # water stress
        cond71 = (cond7 & ((tCCadj < self.Senescence) | (self.tEarlySen > 0)))

        # Early canopy senescence
        cond711 = (cond71 & (self.Ksw['Sen'] < 1))
        self.PrematSenes[cond711] = True
        # No prior early senescence
        cond7111 = (cond711 & (self.tEarlySen == 0))
        self.CCxEarlySen[cond7111] = self.CCprev[cond7111]

        # Increment early senescence GDD counter
        self.tEarlySen[cond711] += dtCC[cond711]

        # Adjust canopy decline coefficient for water stress
        self.water_stress(meteo, beta = False)
        cond7112 = (cond711 & (self.Ksw['Sen'] > 0.99999))
        CDCadj[cond7112] = 0.0001
        cond7113 = (cond711 & np.logical_not(cond7112))
        CDCadj[cond7113] = ((1 - (self.Ksw['Sen'] ** 8)) * self.CDC)[cond7113]

        # Get new canopy cover size after senescence
        cond7114 = (cond711 & (self.CCxEarlySen < 0.001))
        CCsen[cond7114] = 0

        # Get time required to reach CC at end of previous day, given CDCadj
        cond7115 = (cond711 & np.logical_not(cond7114))
        tReq = self.canopy_cover_required_time(self.CCprev, self.CC0adj,
                                               self.CCxEarlySen, self.CGC,
                                               CDCadj, dtCC,
                                               tCCadj, 'CDC')

        # Calculate GDD's for canopy decline and determine new canopy size
        tmp_tCC = tReq + dtCC
        tmp_CCsen = self.canopy_cover_development(self.CC0adj, self.CCxEarlySen, self.CGC,
                                           CDCadj, tmp_tCC, 'Decline')
        CCsen[cond7115] = tmp_CCsen[cond7115]

        
        # Update canopy cover size
        cond7116 = (cond711 & (tCCadj < self.Senescence))

        # Limit CC to CCx
        CCsen[cond7116] = np.clip(CCsen, None, self.CCx)[cond7116]

        # CC cannot be greater than value on previous day
        self.CC[cond7116] = np.clip(self.CC, None, self.CCprev)[cond7116]

        # Update maximum canopy cover size during growing season
        self.CCxAct[cond7116] = self.CC[cond7116]

        # Update CC0 if current CC is less than initial canopy cover size at
        # planting
        cond71161 = (cond7116 & (self.CC < self.CC0))
        self.CC0adj[cond71161] = self.CC[cond71161]
        cond71162 = (cond7116 & np.logical_not(cond71161))
        self.CC0adj[cond71162] = self.CC0[cond71162]

        # Update CC to account for canopy cover senescence due to water stress
        cond7117 = (cond711 & np.logical_not(cond7116))
        self.CC[cond7117] = np.clip(self.CC, None, CCsen)[cond7117]

        # Check for crop growth termination
        cond7118 = (cond711
                    & ((self.CC < 0.001) & np.logical_not(self.CropDead)))
        self.CC[cond7118] = 0
        self.CropDead[cond7118] = True

        # Otherwise there is no water stress
        cond712 = (cond71 & np.logical_not(cond711))
        self.PrematSenes[cond712] = False

        # Rewatering of canopy in late season: get adjusted values of CCx and
        # CDC and update CC
        cond7121 = (cond712 & ((tCCadj > self.Senescence) & (self.tEarlySen > 0)))
        tmp_tCC = tCCadj - dtCC - self.Senescence
        tmp_CCxAdj,tmp_CDCadj = self.update_CCx_and_CDC(self.CCprev, self.CDC, self.CCx, tmp_tCC)
        CCxAdj[cond7121] = tmp_CCxAdj[cond7121]
        CDCadj[cond7121] = tmp_CDCadj[cond7121]
        tmp_tCC = tCCadj - self.Senescence
        tmp_CC = self.canopy_cover_development(self.CC0adj, CCxAdj, self.CGC,
                                        CDCadj, tmp_tCC, 'Decline')
        self.CC[cond7121] = tmp_CC[cond7121]

        # Check for crop growth termination
        cond71211 = (cond7121
                     & ((self.CC < 0.001) & np.logical_not(self.CropDead)))
        self.CC[cond71211] = 0
        self.CropDead[cond71211] = True

        # Reset early senescence counter
        self.tEarlySen[cond712] = 0

        # Adjust CCx for effects of withered canopy
        self.CCxW[cond71] = np.clip(self.CCxW, self.CC, None)[cond71]

        # ######################################################################
        # Calculate canopy size adjusted for micro-advective effects
        
        # Check to ensure potential CC is not slightly lower than actual
        cond8 = (growing_season_index & (self.CC_NS < self.CC))
        self.CC_NS[cond8] = self.CC[cond8]

        cond81 = (cond8 & (tCC < CanopyDevEnd))
        self.CCxAct_NS[cond81] = self.CC_NS[cond81]

        # Actual (with water stress)
        self.CCadj[growing_season_index] = (
            (1.72 * self.CC)
            - (self.CC ** 2)
            + (0.3 * (self.CC ** 3)))[growing_season_index]

        # Potential (without water stress)
        self.CCadj_NS[growing_season_index] = (
            (1.72 * self.CC_NS)
            - (self.CC_NS ** 2)
            + (0.3 * (self.CC_NS ** 3)))[growing_season_index]

        # No canopy outside growing season - set values to zero
        self.CC[np.logical_not(growing_season_index)] = 0
        self.CCadj[np.logical_not(growing_season_index)] = 0
        self.CC_NS[np.logical_not(growing_season_index)] = 0
        self.CCadj_NS[np.logical_not(growing_season_index)] = 0
        self.CCxW[np.logical_not(growing_season_index)] = 0
        self.CCxAct[np.logical_not(growing_season_index)] = 0
        self.CCxW_NS[np.logical_not(growing_season_index)] = 0
        self.CCxAct_NS[np.logical_not(growing_season_index)] = 0
    
    def evap_layer_water_content(self):
        """Function to get water contents in the evaporation layer"""

        # Expand soil properties to compartments
        th_s = self.soil.th_s[self.soil.layerIndex,:]
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]
        th_dry = self.soil.th_dry[self.soil.layerIndex,:]

        # Add rotation, lat, lon dimensions to dz and dzsum
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = self.soil.dz[:,None,None,None] * arr_ones
        dzsum = self.soil.dzsum[:,None,None,None] * arr_ones
        
        # Find compartments covered by evaporation layer
        comp_sto = ((dzsum - dz) < self.EvapZ)
        factor = 1 - ((dzsum - self.EvapZ) / dz)
        factor = np.clip(factor, 0, 1) * growing_season_index * comp_sto
            
        # Water storages in evaporation layer (mm)
        self.Wevap['Act'] = np.sum((factor * 1000 * self.th * dz), axis=0)
        self.Wevap['Act'] = np.clip(self.Wevap['Act'], 0, None)
        self.Wevap['Sat'] = np.sum((factor * 1000 * th_s * dz), axis=0)
        self.Wevap['Fc'] = np.sum((factor * 1000 * th_fc * dz), axis=0)
        self.Wevap['Wp'] = np.sum((factor * 1000 * th_wp * dz), axis=0)
        self.Wevap['Dry'] = np.sum((factor * 1000 * th_dry * dz), axis=0)
                    
    def soil_evaporation(self, rotation, meteo, currTimeStep):
        """Function to calculate daily soil evaporation in AOS"""

        # Add dimensions to dz,dzsum
        arr_ones = np.ones((nr, nlat, nlon))[None,:,:,:]
        dz = self.soil.dz[:,None,None,None] * arr_ones
        dzsum = self.soil.dzsum[:,None,None,None] * arr_ones

        # Add rotation dimension to meteo vars
        et0 = meteo.referencePotET[None,:,:] * np.ones((self.nRotation))[:,None,None]
        prec = meteo.precipitation[None,:,:] * np.ones((self.nRotation))[:,None,None]

        # ######################################################################
        # Prepare stage 2 evaporation (REW gone) (Only do this if it is first
        # day of the simulation, or if it is first day of growing season and
        # not simulating off-season
        cond1 = (
            (currTimeStep.timeStepPCR == 1)
            | ((self.DAP == 1) & (self.OffSeason == False)))
                 
        # Reset storage in surface soil layer to zero
        self.Wsurf[cond1] = 0
        # Set evaporation depth to minimum
        self.EvapZ[cond1] = self.soil.EvapZmin[cond1]
        # Trigger stage 2 evaporation
        self.Stage2[cond1] = True
        # Get relative water content for start of stage 2 evaporation
        Wevap = self.evap_layer_water_content(Wevap, MASK = cond1)
        self.Wstage2[cond1] = (
            (Wevap['Act'] - (Wevap['Fc'] - self.soil.REW))
            / (Wevap['Sat'] - (Wevap['Fc'] - self.soil.REW)))[cond1]
        self.Wstage2[cond1] = (np.round((self.Wstage2 * 100)) / 100)[cond1]
        self.Wstage2[cond1] = np.clip(self.Wstage2, 0, None)[cond1]


        # ######################################################################
        # Prepare soil evaporation stage 1: adjust water in surface evaporation
        # layer for any infiltration (only do this if rainfall occurs or when
        # irrigation is triggered)
        cond2 = ((prec > 0) | np.any(((self.Irr > 0) & (self.rotation.IrrMethod != 4)), axis=0))

        # Update storage in surface evaporation layer for incoming infiltration
        cond21 = (cond2 & (self.Infl > 0))
        self.Wsurf[cond21] = self.Infl[cond21]

        # Water stored in surface evaporation layer cannot exceed REW
        self.Wsurf[cond21] = np.clip(self.Wsurf, None, self.soil.REW)[cond21]

        # Reset variables
        self.Wstage2[cond21] = 0  # TODO: is this right?
        self.EvapZ[cond21] = self.soil.EvapZmin[cond21]
        self.Stage2[cond21] = False

        # ######################################################################
        # Calculate potential soil evaporation rate (mm/day)

        # NB calculation of potential evaporation rate uses crop-specific
        # parameters
        
        # Adjust time for any delayed development
        if rotation.CalendarType == 1:
            tAdj = (rotation.DAP - rotation.DelayedCDs) * growing_season_index
        elif rotation.CalendarType == 2:
            tAdj = (rotation.GDDcum - rotation.DelayedGDDs) * growing_season_index
            
        # Calculate maximum potential soil evaporation and potential soil
        # evaporation given current canopy size
        EsPotMax = (self.soil.Kex * et0 * (1 - self.CCxW * (self.soil.fwcc / 100))) * growing_season_index
        EsPot = (self.soil.Kex * (1 - self.CCadj) * et0) * growing_season_index

        # Adjust potential soil evaporation for effects of withered canopy
        cond3 = (growing_season_index & (tAdj > self.rotation.Senescence) & (self.CCxAct > 0))
        mult = np.zeros((nr, nlat, nlon))
        cond31 = (cond3 & (self.CC > (self.CCxAct / 2)))
        cond311 = (cond31 & (self.CC <= self.CCxAct))
        mult[cond311] = ((self.CCxAct - self.CC) / (self.CCxAct / 2))[cond311]
        cond32 = (cond3 & np.logical_not(cond31))
        mult[cond32] = 1
        EsPot[cond3] = (EsPot * (1 - self.CCxAct * (self.soil.fwcc / 100) * mult))[cond3]
        CCxActAdj = ((1.72 * self.CCxAct) + (self.CCxAct ** 2) - 0.3 * (self.CCxAct ** 3))
        EsPotMin = np.zeros((nr, nlat, nlon))
        EsPotMin[cond3] = (self.soil.Kex * (1 - CCxActAdj) * et0)[cond3]
        EsPotMin = np.clip(EsPotMin, 0, None)
        EsPot[cond3] = np.clip(EsPot, EsPotMin, EsPotMax)[cond3]

        cond4 = (growing_season_index & self.PrematSenes)
        EsPot[cond4] = np.clip(EsPot, None, EsPotMax)[cond4]

        # No canopy cover outside of growing season so potential soil
        # evaporation only depends on reference evapotranspiration
        cond5 = np.logical_not(growing_season_index)
        EsPot[cond5] = (self.soil.Kex * et0)[cond5]

        # ######################################################################
        # Adjust potential soil evaporation for mulches and/or partial wetting
        
        # self.rotation.Mulches
        EsPotMul = np.zeros((nr, nlat, nlon))
        cond5 = (self.SurfaceStorage < 0.000001)

        # No mulches present
        cond51 = (cond5 & (self.rotation.Mulches == 0))
        EsPotMul[cond51] = EsPot[cond51]

        # self.rotation.Mulches present (percentage soil surface covered may vary depending
        # on whether within or outside growing season"
        cond52 = (cond5 & (self.rotation.Mulches == 1))
        cond521 = (cond52 & growing_season_index)
        EsPotMul[cond521] = (EsPot * (1 - self.rotation.fMulch * (self.rotation.MulchPctGS / 100)))[cond521]
        cond522 = (cond52 & np.logical_not(growing_season_index))
        EsPotMul[cond522] = (EsPot * (1 - self.rotation.fMulch * (self.rotation.MulchPctOS / 100)))[cond522]

        # Surface is flooded - no adjustment of potential soil evaporation for
        # mulches
        cond6 = np.logical_not(cond5)
        EsPotMul[cond6] = EsPot[cond6]

        # ######################################################################
        # Partial surface wetting by irrigation

        # Only apply adjustment if irrigation occurs and not in net irrigation
        # mode 
        EsPotIrr = np.zeros((nr, nlat, nlon))
        cond7 = ((self.Irr > 0) & (self.rotation.IrrMethod != 4))
        cond71 = (cond7 & ((prec > 0) | (self.SurfaceStorage > 0)))
        EsPotIrr[cond71] = EsPot[cond71]
        cond72 = (cond7 & np.logical_not(cond71))
        EsPotIrr[cond72] = (EsPot * (self.rotation.WetSurf / 100))[cond72]

        # Otherwise do not adjust for partial surface wetting
        cond8 = np.logical_not(cond7)
        EsPotIrr[cond8] = EsPot[cond8]

        # Assign minimum value (mulches and partial wetting don't combine)
        EsPot = np.minimum(EsPotIrr, EsPotMul)

        # ######################################################################
        # Surface evaporation
        
        # Initialise actual evaporation counter
        self.EsAct = np.zeros((nr, nlat, nlon))

        # Evaporate surface storage
        cond9 = (self.SurfaceStorage > 0)
        cond91 = (cond9 & (self.SurfaceStorage > EsPot))
        self.EsAct[cond91] = EsPot[cond91]
        self.SurfaceStorage[cond91] = (self.SurfaceStorage - self.EsAct)[cond91]

        # Otherwise surface storage is not sufficient to meet all potential
        # soil evaporation
        cond92 = (cond9 & np.logical_not(cond91))
        self.EsAct[cond92] = self.SurfaceStorage[cond92]

        # Update surface storage, evaporation layer depth, stage
        self.SurfaceStorage[cond92] = 0
        self.Wsurf[cond92] = self.soil.REW[cond92]
        self.Wstage2[cond92] = 0 # not sure if this is correct
        self.EvapZ[cond92] = self.soil.EvapZmin[cond92]
        self.Stage2[cond92] = False

        # ######################################################################
        # Stage 1 evaporation

        # Determine total water to be extracted
        ToExtract = EsPot - self.EsAct

        # Determine total water to be extracted in stage one (limited by
        # surface layer water storage)
        ExtractPotStg1 = np.minimum(ToExtract,self.Wsurf)

        # Extract water
        cond10 = (ExtractPotStg1 > 0)
        
        # Determine fraction of compartments covered by evaporation layer
        comp_sto = ((dzsum - dz) < self.soil.EvapZmin)
        factor = 1 - ((dzsum - self.soil.EvapZmin) / dz)
        factor = np.clip(factor, 0, 1) * comp_sto
        
        comp_sto = np.sum(comp_sto, axis=0)
        comp = 0
        while np.any((comp < comp_sto) & (ExtractPotStg1 > 0)):
            
            cond101 = ((comp < comp_sto) & (ExtractPotStg1 > 0))
            
            # Water available in compartment for extraction (mm)
            Wdry = 1000 * th_dry[comp,:] * dz[comp,:]  
            W = 1000 * self.th[comp,:] * dz[comp,:]
            AvW = np.zeros((nr, nlat, nlon))
            AvW[cond101] = ((W - Wdry) * factor[comp,:])[cond101]
            AvW = np.clip(AvW, 0, None)

            # Determine amount by which to adjust variables
            chng = np.zeros((nr, nlat, nlon))
            cond1011 = (cond101 & (AvW >= ExtractPotStg1))
            chng[cond1011] = ExtractPotStg1[cond1011]
            cond1012 = (cond101 & np.logical_not(cond1011))
            chng[cond1012] = AvW[cond1012]

            # NB no need for index because chng = 0 if cond1011|1012 not met
            self.EsAct += chng       # actual evaporation 
            W -= chng           # depth of water in current compartment
            ToExtract -= chng   # total water to be extracted
            ExtractPotStg1 -= chng  # water to be extracted from surface layer

            # Update water content
            self.th[comp,:][cond101] = (W / (1000 * dz[comp,:]))[cond101]
            comp += 1

        # Update surface evaporation layer water balance
        self.Wsurf[cond10] -= self.EsAct[cond10]
        cond102 = (cond10 & ((self.Wsurf < 0) | (ExtractPotStg1 > 0.0001)))
        self.Wsurf[cond102] = 0

        # If surface storage completely depleted, prepare stage 2 evaporation
        cond103 = (cond10 & (self.Wsurf < 0.0001))
        # Get water contents
        Wevap = self.evap_layer_water_content(Wevap, MASK = cond103)

        # Proportional water storage for start of stage two evaporation
        self.Wstage2[cond103] = (
            (Wevap['Act'] - (Wevap['Fc'] - self.soil.REW)) /
            (Wevap['Sat'] - (Wevap['Fc'] - self.soil.REW)))[cond103]
        self.Wstage2[cond103] = (np.round((self.Wstage2 * 100)) / 100)[cond103]
        self.Wstage2[cond103] = np.clip(self.Wstage2, 0, None)[cond103]

        # ######################################################################
        # Stage 2 evaporation

        # Extract water
        cond11 = (ToExtract > 0)
        if np.any(cond11):

            # Start stage 2
            self.Stage2[cond11] = True

            # Get sub-daily evaporative demand
            self.EvapTimeSteps = 20  # TODO: add to options
            Edt = ToExtract / self.EvapTimeSteps

            # Loop sub-daily time steps
            for jj in range(self.EvapTimeSteps):

                # Get current water storage
                Wevap = self.evap_layer_water_content(Wevap, MASK = cond11)
                
                # Get water storage (mm) at start of stage 2 evaporation
                Wupper = (
                    self.Wstage2
                    * (Wevap['Sat'] - (Wevap['Fc'] - self.soil.REW))
                    + (Wevap['Fc'] - self.soil.REW))
                
                # Get water storage (mm) when there is no evaporation
                Wlower = Wevap['Dry']
                
                # Get relative depletion of evaporation storage in stage 2
                Wrel = ((Wevap['Act'] - Wlower) / (Wupper - Wlower))
                
                # Check if need to expand evaporative layer
                cond111 = (cond11 & (self.soil.EvapZmax > self.soil.EvapZmin))
                Wcheck = (
                    self.soil.fWrelExp
                    * ((self.soil.EvapZmax - self.EvapZ)
                       / (self.soil.EvapZmax - self.soil.EvapZmin)))
                
                while np.any(cond111 & (Wrel < Wcheck) & (self.EvapZ < self.soil.EvapZmax)):
                    cond1111 = (cond111 & (Wrel < Wcheck) & (self.EvapZ < self.soil.EvapZmax))

                    # Expand evaporation layer by 1mm
                    self.EvapZ[cond111] += 0.001

                    # Recalculate current water storage for new EvapZ
                    Wevap = evap_layer_water_content(Wevap)
                    Wupper = (
                        self.Wstage2
                        * (Wevap['Sat'] - (Wevap['Fc'] - self.soil.REW))
                        + (Wevap['Fc'] - self.soil.REW))
                    Wlower = Wevap['Dry']
                    Wrel = ((Wevap['Act'] - Wlower) / (Wupper - Wlower))
                    Wcheck = (
                        self.soil.fWrelExp
                        * ((self.soil.EvapZmax - self.EvapZ)
                           / (self.soil.EvapZmax - self.soil.EvapZmin)))

                # Get stage 2 evaporation reduction coefficient
                Kr = ((np.exp(self.soil.fevap * Wrel) - 1) / (np.exp(self.soil.fevap) - 1))
                Kr = np.clip(Kr, None, 1)

                # Get water to extract (NB Edt is zero in cells which do not
                # need stage 2, so no need for index)
                ToExtractStg2 = (Kr * Edt)

                # Determine fraction of compartments covered by evaporation layer
                comp_sto = ((dzsum - dz) < self.soil.EvapZmin)
                factor = 1 - ((dzsum - self.soil.EvapZmin) / dz)
                factor = np.clip(factor, 0, 1) * comp_sto

                comp_sto = np.sum(comp_sto, axis=0)
                comp = 0
                while np.any(cond11 & (comp < comp_sto) & (ToExtractStg2 > 0)):

                    cond111 = (cond11 & (comp < comp_sto) & (ToExtractStg2 > 0))
                    
                    # Water available in compartment for extraction (mm)
                    Wdry = 1000 * th_dry[comp,:] * dz[comp,:]  
                    W = 1000 * self.th[comp,:] * dz[comp,:]
                    AvW = np.zeros((nr, nlat, nlon))
                    AvW[cond111] = ((W - Wdry) * factor[comp,:])[cond111]
                    AvW = np.clip(AvW, 0, None)

                    # Determine amount by which to adjust variables
                    chng = np.zeros((nr, nlat, nlon))
                    cond1111 = (cond111 & (AvW >= ExtractPotStg1))
                    chng[cond1111] = ExtractPotStg1[cond1111]
                    cond1112 = (cond111 & np.logical_not(cond1111))
                    chng[cond1112] = AvW[cond1112]

                    # NB no need for index because chng = 0 if cond1111|1112 not met
                    self.EsAct += chng       # actual evaporation 
                    W -= chng           # depth of water in current compartment
                    ToExtract -= chng   # total water to be extracted
                    ExtractPotStg1 -= chng  # water to be extracted from surface layer

                    # Update water content
                    self.th[comp,:][cond111] = (W / (1000 * dz[comp,:]))[cond111]
                    comp += 1
                
        # ######################################################################
        # Store potential evaporation for irrigation calculations on next day

        self.Epot = EsPot

    def aeration_stress(self):
        """Function to calculate aeration stress coefficient"""

        # Determine aeration stress (root zone)
        cond1 = (self.thRZ['Act'] > self.thRZ['Aer'])

        # Calculate aeration stress coefficient
        self.Ksa_Aer = np.ones((self.nRotation, self.nLat, self.nLon))
        cond11 = (cond1 & (self.AerDays < self.LagAer))
        stress = (1 - ((self.thRZ['Sat'] - self.thRZ['Act'])
                       / (self.thRZ['Sat'] - self.thRZ['Aer'])))
        self.Ksa_Aer[cond11] = (1 - ((self.AerDays / 3) * stress))[cond11]
        cond12 = (cond1 & np.logical_not(cond11))
        self.Ksa_Aer[cond12] = ((self.thRZ['Sat'] - self.thRZ['Act'])
                                / (self.thRZ['Sat'] - self.thRZ['Aer']))[cond12]

        # Increment aeration days counter, or set to zero if there is no stress
        self.AerDays[cond1] += 1
        self.AerDays[np.logical_not(cond1)] = 0
        self.AerDays = np.clip(self.AerDays, None, self.LagAer)
        
    def transpiration(self, meteo):
        """Function to calculate crop transpiration on current day"""

        # Expand soil properties to compartments
        th_s = self.soil.th_s[self.soil.layerIndex,:]
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]
        th_dry = self.soil.th_dry[self.soil.layerIndex,:]

        # Add dimensions to dz,dzsum
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = self.soil.dz[:,None,None,None] * arr_ones
        dzsum = self.soil.dzsum[:,None,None,None] * arr_ones
        
        # Add rotation dimension to ET0
        et0 = meteo.referencePotET[None,:,:] * np.ones((self.nRotation))[:,None,None]

        # ######################################################################
        # Calculate transpiration (if in growing season)

        # 1. No prior water stress

        # Update ageing days counter
        cond1 = (growing_season_index & (self.DAP > self.MaxCanopyCD))
        self.AgeDays_NS[cond1] = (self.DAP - self.MaxCanopyCD)[cond1]

        # Update crop coefficient for ageing of canopy
        Kcb_NS = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond2 = (growing_season_index & (self.AgeDays_NS > 5))
        Kcb_NS[cond2] = (self.Kcb - ((self.AgeDays_NS - 5) * (self.fage / 100))
                         * self.CCxW_NS)[cond2]
        cond3 = (growing_season_index & np.logical_not(cond2))
        Kcb_NS[cond3] = self.Kcb[cond3]

        # Update crop coefficient for C02 concentration ***TODO***
        # cond4 = (growing_season_index & (CurrentConc > RefConc))
        # Kcb_NS[cond4] *= (1 - 0.05 * (CurrentConc - RefConc) / (550 - RefConc))[cond4]

        # Determine potential transpiration rate (no water stress)
        self.TrPot_NS = Kcb_NS * self.CCadj_NS * et0 * growing_season_index

        # 2. Potential prior water stress and/or delayed development

        # Update ageing days counter
        DAPadj = (self.DAP - self.DelayedCDs) * growing_season_index
        cond5 = (growing_season_index & (DAPadj > self.MaxCanopyCD))
        self.AgeDays[cond5] = (DAPadj - self.MaxCanopyCD)[cond5]

        # Update crop coefficient for ageing of canopy
        Kcb = self.Kcb
        cond6 = (growing_season_index & (self.AgeDays > 5))
        Kcb[cond6] = (self.Kcb - ((self.AgeDays - 5) * (self.fage / 100)) * self.CCxW)[cond6]

        # Update crop coefficient for CO2 concentration ***TODO***
        # cond8 = (growing_season_index & (CurrentConc > RefConc))
        # Kcb[cond8] *= (1 - 0.05 * (CurrentConc - RefConc) / (550 - RefConc))[cond4][cond8]

        # Determine potential transpiration rate
        self.TrPot0 = Kcb * (self.CCadj) * et0 * growing_season_index

        # Correct potential transpiration for dying green canopy effects
        cond9 = (growing_season_index & (self.CC < self.CCxW))
        cond91 = (cond9 & (self.CCxW > 0.001) & (self.CC > 0.001))
        self.TrPot0[cond91] *= ((self.CC / self.CCxW) ** self.a_Tr)[cond91]

        # ######################################################################
        # Calculate surface layer transpiration

        cond10 = (growing_season_index
                  & (self.SurfaceStorage > 0)
                  & (self.DaySubmerged < self.LagAer))

        # Initialise variables
        TrPot = np.zeros((self.nRotation, self.nLat, self.nLon))
        
        # Update submergence days counter
        self.DaySubmerged[cond10] += 1

        # Update anaerobic conditions counter for each compartment
        cond10_comp = np.broadcast_to(cond10, self.AerDaysComp.shape)
        self.AerDaysComp[cond10_comp] += 1 
        LagAer_comp = np.broadcast_to(self.LagAer, self.AerDaysComp.shape)
        self.AerDaysComp[cond10_comp] = np.clip(self.AerDaysComp, None, LagAer_comp)[cond10_comp]

        # Reduce actual transpiration that is possible to account for aeration
        # stress due to extended submergence
        fSub = 1 - (self.DaySubmerged / LagAer)

        # Transpiration occurs from surface storage
        cond101 = (cond10 & (self.SurfaceStorage > (fSub * self.TrPot0)))
        # self.SurfaceStorage[cond101] -= (fSub * self.TrPot0)[cond101]
        self.TrAct0[cond101] = (fSub * self.TrPot0)[cond101]

        # Otherwise there is no transpiration from surface storage
        cond102 = (cond10 & np.logical_not(cond101))
        self.TrAct0[cond102] = 0

        # More water can be extracted from soil profile for transpiration
        cond103 = (cond10 & (self.TrAct0 < (fSub * self.TrPot0)))
        TrPot[cond103] = ((fSub * self.TrPot0) - self.TrAct0)[cond103]

        # ######################################################################
        # "Update potential root zone transpiration for water stress"

        # "Determine root zone water content"
        self.root_zone_water()

        # Calculate water stress coefficients
        self.water_stress(meteo, beta = True)

        # Calculate aeration stress coefficients
        self.aeration_stress()

        # Maximum stress effect
        Ks = np.minimum(self.Ksw['StoLin'], self.Ksa_Aer)

        # Update potential transpiration in root zone
        cond11 = (growing_season_index & (self.IrrMethod != 4))
        TrPot[cond11] *= Ks[cond11]

        # ######################################################################
        # Determine compartments covered by root zone

        rootdepth = np.maximum(self.Zmin, self.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        comp_sto = ((dzsum - dz) < rootdepth)
        
        # Fraction of compartment covered by root zone (zero in compartments
        # NOT covered by the root zone)
        RootFact = 1 - ((dzsum - rootdepth) / dz)
        RootFact = np.clip(RootFact, 0, 1) * comp_sto

        # ######################################################################
        # Determine maximum sink term for each compartment

        SxComp = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))

        # Net irrigation mode
        cond12 = (growing_season_index & (self.IrrMethod == 4))
        cond12_comp = np.broadcast_to(cond12, SxComp.shape)
        SxComp[cond12_comp] = np.broadcast_to(((self.SxTop + self.SxBot) / 2), SxComp.shape)[cond12_comp]

        # Otherwise sink term declines linearly with depth
        cond13 = (growing_season_index & np.logical_not(cond12))

        if (np.any(cond13)):
            comp = 0
            comp_sto_sum = np.sum(comp_sto, axis=0)
            SxCompBot = self.SxTop
            while np.any(comp < comp_sto_sum):
                SxCompTop = SxCompBot
                cond131 = (cond13 & (dzsum[comp,:] <= rootdepth))
                SxCompBot[cond131] = (
                    self.SxBot * self.rCor
                    + ((self.SxTop - self.SxBot * self.rCor)
                       * ((rootdepth - dzsum[comp,:]) / rootdepth)))[cond131]
                
                cond132 = (cond13 & np.logical_not(cond131))
                SxCompBot[cond132] = (self.SxBot * self.rCor)[cond132]
                SxComp[cond13] = ((SxCompTop + SxCompBot) / 2)[cond13]
                comp += 1
                
        SxComp *= comp_sto

        # ######################################################################
        # Extract water

        ToExtract = TrPot
        self.TrAct = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond14_ini = (growing_season_index & (ToExtract > 0))
        
        if (np.any(cond14_ini)):
            comp = 0
            comp_sto_sum = np.sum(comp_sto, axis=0)
            while np.any((comp < comp_sto_sum) & (ToExtract > 0)):

                cond14 = (growing_season_index & (comp_sto[comp,:]) & (ToExtract > 0))

                # Determine TAW for compartment
                thTAW = th_fc[comp,:] - th_wp[comp,:]
                p_up_sto = np.ones((self.nRotation, self.nLat, self.nLon))
                cond141 = (cond14 & (ETadj == 1))
                p_up_sto[cond141] = (
                    self.p_up2
                    + (0.04 * (5 - et0))
                    * (np.log10(10 - 9 * self.p_up2)))[cond141]
                # Determine critical water content at which stomatal closure
                # will occur in compartment
                
                thCrit = (th_fc[comp,:] - (thTAW * p_up_sto))

                # Check for soil water stress
                KsComp = np.zeros((self.nRotation, self.nLat, self.nLon))
                # No stress
                cond142 = (cond14 & (self.th[comp,:] >= thCrit))
                KsComp[cond142] = 1
                # Transpiration from compartment is affected by water stress
                cond143 = (cond14 & (self.th[comp,:] > th_wp[comp,:]) & np.logical_not(cond142))
                Wrel = ((th_fc[comp,:] - self.th[comp,:]) / (th_fc[comp,:] - th_wp[comp,:]))
                pRel = ((Wrel - self.p_up2) / (self.p_lo2 - self.p_up2))
                KsComp[cond143] = (1 - ((np.exp(pRel * self.fshape_w2) - 1) / (np.exp(self.fshape_w2) - 1)))[cond143]
                KsComp = np.clip(KsComp, 0, 1)
                KsComp[pRel <= 0] = 1
                KsComp[pRel >= 1] = 0
                
                # Otherwise no transpiration is possible from the compartment
                # as water does not exceed wilting point (no need to explicitly
                # assign zero because KsComp is initialized with zeros)

                # Adjust compartment stress factor for aeration stress
                AerComp = np.zeros((self.nRotation, self.nLat, self.nLon))
                # Full aeration stress - no transpiration possible from
                # compartment
                cond144 = (cond14 & (self.DaySubmerged >= self.LagAer))
                cond145 = (cond14
                           & (self.th[comp,:] > (th_s[comp,:] - (self.Aer / 100)))
                           & np.logical_not(cond144))
                self.AerDaysComp[comp,:][cond145] += 1
                fAer = np.ones((self.nRotation, self.nLat, self.nLon))
                cond1451 = (cond145 & self.AerDaysComp[comp,:] >= self.LagAer)
                self.AerDaysComp[cond1451] = self.LagAer[cond1451]
                fAer[cond1451] = 0

                # Calculate aeration stress factor
                AerComp[cond145] = (
                    (th_s[comp,:] - self.th[comp,:])
                    / (th_s[comp,:] - (th_s[comp,:] - (self.Aer / 100))))[cond145]
                AerComp = np.clip(AerComp, 0, None)
                AerComp[cond145] = (
                    (fAer + (self.AerDaysComp[comp,:] - 1) * AerComp)
                    / (fAer + self.AerDaysComp[comp,:] - 1))[cond145]

                # Otherwise there is no aeration stress as number of submerged
                # days does not exceed threshold for initiation of aeration
                # stress
                cond146 = (cond14 & np.logical_not(cond144 | cond145))
                AerComp[cond146] = 1
                self.AerDaysComp[comp,:][cond146] = 0

                # Extract water
                ThToExtract = ((ToExtract / 1000) / dz[comp,:])
                Sink = np.zeros((self.nRotation, self.nLat, self.nLon))
                # Don't reduce compartment sink for stomatal water stress if in
                # net irrigation mode. Stress only occurs due to deficient
                # aeration conditions
                cond147 = (cond14 & self.IrrMethod == 4)
                Sink[cond147] = (AerComp * SxComp[comp,:] * RootFact[comp,:])[cond147]

                # Otherwise, reduce compartment sink for greatest of stomatal
                # and aeration stress
                cond148 = (cond14 & np.logical_not(cond147))
                cond1481 = (cond148 & (KsComp == AerComp))
                Sink[cond1481] = (KsComp * SxComp[comp,:] * RootFact[comp,:])[cond1481]
                cond1482 = (cond148 & np.logical_not(cond1481))
                Sink[cond1482] = (np.minimum(KsComp,AerComp) * SxComp[comp,:] * RootFact[comp,:])[cond1482]

                # Limit extraction to demand
                Sink = np.clip(Sink, None, ThToExtract)

                # Limit extraction to avoid compartment water content dropping
                # below air dry
                cond149 = (cond14 & ((self.th[comp,:] - Sink) < th_dry[comp,:]))
                Sink[cond149] = (self.th[comp,:] - th_dry[comp,:])[cond149]
                Sink = np.clip(Sink, 0, None)

                # Update water content in compartment
                self.th[comp,:][cond14] -= Sink[cond14]

                # Update amount of water to extract
                ToExtract[cond14] -= (Sink * 1000 * dz[comp,:])[cond14]

                # Update actual transpiration
                self.TrAct[cond14] += (Sink * 1000 * dz[comp,:])[cond14]

                # Update compartment counter
                comp += 1
            
        # ######################################################################
        # Add net irrigation water requirement (if this mode is specified)
        cond15 = (growing_season_index & (self.IrrMethod == 4) & (TrPot > 0))
        self.IrrNet = np.zeros((self.nRotation, self.nLat, self.nLon))
        self.root_zone_water()
        thCrit = self.thRZ['Wp'] + ((self.NetIrrSMT / 100) * (self.thRZ['Fc'] - self.thRZ['Wp']))
        cond151 = (cond15 & (self.thRZ['Act'] < thCrit))
        cond151_comp = np.broadcast_to(cond151, self.th.shape)

        # Calculate thCrit in each compartment
        thCrit_comp = (th_wp + ((self.NetIrrSMT / 100) * (th_fc - th_wp)))
        # Determine necessary change in water content in compartments to reach
        # critical water content
        dWC = RootFact * (thCrit_comp - self.th * 1000 * dz)
        self.th[cond151_comp] = (self.th + (dWC / (1000 * dz)))[cond151_comp]
        self.IrrNet[cond151] = np.sum(dWC, axis=0)[cond151]

        # Update net irrigation counter for the growing season
        self.IrrNetCum += self.IrrNet

        # # No net irrigation as potential transpiration is zero
        # cond16 = (growing_season_index & (self.IrrMethod == 4) & (TrPot <= 0))
        # IrrNet[cond16] = 0

        # # No net irrigation as not in net irrigation mode
        # cond17 = (growing_season_index & np.logical_not(cond15 | cond16))
        # IrrNet[cond17] = 0
        # self.IrrNetCum[cond17] = 0

        # ######################################################################
        # Add any surface transpiration to root zone total
        self.TrAct += self.TrAct0

        # ######################################################################
        # Feedback with canopy cover development

        # If actual transpiration is zero then no canopy cover growth can occur
        cond16 = (growing_season_index & ((self.CC - self.CCprev) > 0.005) & (self.TrAct > 0))
        self.CC[cond16] = self.CCprev[cond16]

        # ######################################################################
        # Update transpiration ratio
        cond17 = (growing_season_index & (self.TrPot0 > 0))
        cond171 = (cond17 & (self.TrAct < self.TrPot0))
        self.TrRatio[cond171] = (self.TrAct / self.TrPot0)[cond171]
        cond172 = (cond17 & np.logical_not(cond171))
        self.TrRatio[cond172] = 1
        cond18 = (growing_season_index & np.logical_not(cond17))
        self.TrRatio[cond18] = 1
        self.TrRatio = np.clip(self.TrRatio, 0, 1)

        # No transpiration or irrigation if outside growing season
        self.TrAct[np.logical_not(growing_season_index)] = 0
        self.TrPot0[np.logical_not(growing_season_index)] = 0
        self.TrPot_NS[np.logical_not(growing_season_index)] = 0
        self.IrrNet[np.logical_not(growing_season_index)] = 0
        self.IrrNetCum[np.logical_not(growing_season_index)] = 0

        # Store potential transpiration for irrigation calculations on next day
        self.Tpot = self.TrPot0

    def groundwater_inflow(self, groundwater):
        """Function to calculate capillary rise in the presence of a 
        shallow groundwater table
        """
        # Expand soil properties to compartments
        th_s = self.soil.th_s[self.soil.layerIndex,:]
        
        # Add dimensions to dz,dzsum
        arr_ones = np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]
        dz = self.soil.dz[:,None,None,None] * arr_ones
        dzsum = self.soil.dzsum[:,None,None,None] * arr_ones

        # Initialize groudwater inflow array
        self.GwIn = np.zeros((self.nRotation, self.nLat, self.nLon))
        
        # Water table in soil profile: calculate horizontal inflow; get
        # groundwater table elevation on current day
        zBot = np.cumsum(dz, axis=0)
        zTop = zBot - dz
        zMid = (zTop + zBot) / 2

        # For compartments below water table, set to saturation
        dth = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))
        cond1 = (self.WTinSoil & (zMid >= self.zGW))
        cond11 = (cond1 & (self.th < th_s))

        # Update water content
        dth[cond11] = (th_s - self.th)[cond11]
        self.th[cond11] = th_s[cond11]

        # Update groundwater inflow
        GwIn_comp = dth * 1000 * dz
        self.GwIn = np.sum(GwIn_comp, axis=0)

    def harvest_index_ref_current_day(self):
        """Function to calculate reference (no adjustment for stress 
        effects) harvest index on current day
        """
        # Check if in yield formation period
        if self.CalendarType == 1:
            tAdj = self.DAP - self.DelayedCDs
        elif self.CalendarType == 2:
            tAdj = self.GDDcum - self.DelayedGDDs

        self.YieldForm = (growing_season_index & (tAdj > self.HIstart))

        # Get time for harvest index calculation
        HIt = self.DAP - self.DelayedCDs - self.HIstartCD - 1

        # Yet to reach time for HI build-up
        cond1 = (growing_season_index & (HIt <= 0))
        self.HIref[cond1] = 0
        self.PctLagPhase[cond1] = 0

        cond2 = (growing_season_index & np.logical_not(cond1))

        # HI cannot develop further as canopy is too small (no need to do
        # anything here as HIref doesn't change)
        cond21 = (cond2 & (self.CCprev <= (self.CCmin * self.CCx)))
        cond22 = (cond2 & np.logical_not(cond21))

        # If crop type is leafy vegetable or root/tuber then proceed with
        # logistic growth (i.e. no linear switch)
        cond221 = (cond22 & ((self.CropType == 1) | (self.CropType == 2)))
        self.PctLagPhase[cond221] = 100
        self.HIref[cond221] = ((self.HIini * self.HI0) / (self.HIini + (self.HI0 - self.HIini) * np.exp(-self.HIGC * HIt)))[cond221]
        # Harvest index approaching maximum limit
        cond2211 = (cond221 & (self.HIref >= (0.9799 * HI0)))
        self.HIref[cond2211] = HI0[cond2211]

        cond222 = (cond22 & (self.CropType == 3))
        # Not yet reached linear switch point, therefore proceed with logistic
        # build-up
        cond2221 = (cond222 & (HIt < self.tLinSwitch))
        self.PctLagPhase[cond2221] = (100 * (HIt / self.tLinSwitch))[cond2221]
        self.HIref[cond2221] = ((self.HIini * self.HI0) / (self.HIini + (self.HI0 - self.HIini) * np.exp(-self.HIGC * HIt)))[cond2221]
        cond2222 = (cond222 & np.logical_not(cond2221))
        self.PctLagPhase[cond2222] = 100
        # Calculate reference harvest index for current day (logistic portion)
        self.HIref[cond2222] = ((self.HIini * self.HI0) / (self.HIini + (self.HI0 - self.HIini) * np.exp(-self.HIGC * self.tLinSwitch)))[cond2222]
        # Calculate reference harvest index for current day (total = logistic + linear)
        self.HIref[cond2222] += (self.dHILinear * (HIt - self.tLinSwitch))[cond222]

        # Limit HIref and round off computed value
        cond223 = (cond22 & (self.HIref > HI0))
        self.HIref[cond223] = self.HI0[cond223]
        cond224 = (cond22 & (self.HIref <= (self.HIini + 0.004)))
        self.HIref[cond224] = 0
        cond225 = (cond22 & ((self.HI0 - self.HIref) < 0.004))
        self.HIref[cond225] = self.HI0[cond225]

        # Reference harvest index is zero outside growing season
        self.HIref[np.logical_not(growing_season_index)] = 0

    def temperature_stress(self, meteo):
        """Function to calculate temperature stress coefficients"""

        # Add rotation dimension to meteo vars
        tmin = meteo.tmin[None,:,:] * np.ones((self.nRotation))[:,None,None]
        tmax = meteo.tmax[None,:,:] * np.ones((self.nRotation))[:,None,None]
        
        # Calculate temperature stress coefficient affecting biomass growth
        KsBio_up = 1
        KsBio_lo = 0.02
        fshapeb = -1 * (np.log(((KsBio_lo * KsBio_up) - 0.98 * KsBio_lo) / (0.98 * (KsBio_up - KsBio_lo))))

        # Calculate temperature stress effects on biomass production
        cond1 = (self.BioTempStress == 0)
        self.Kst['Bio'][cond1] = 1
        cond2 = (self.BioTempStress == 1)
        cond21 = (cond2 & (self.GDD >= self.GDD_up))
        self.Kst['Bio'][cond21] = 1
        cond22 = (cond2 & (self.GDD <= self.GDD_lo))
        self.Kst['Bio'][cond22] = 0
        cond23 = (cond2 & np.logical_not(cond21 | cond22))
        GDDrel = (self.GDD - self.GDD_lo) / (self.GDD_up - self.GDD_lo)
        self.Kst['Bio'][cond23] = ((KsBio_up * KsBio_lo) / (KsBio_lo + (KsBio_up - KsBio_lo) * np.exp(-fshapeb * GDDrel)))[cond23]
        self.Kst['Bio'][cond23] = (self.Kst['Bio'] - KsBio_lo * (1 - GDDrel))[cond23]

        # Calculate temperature stress coefficients affecting crop pollination
        KsPol_up = 1
        KsPol_lo = 0.001

        # Calculate effects of heat stress on pollination
        cond3 = (self.PolHeatStress == 0)
        self.Kst['PolH'][cond3] = 1
        cond4 = (self.PolHeatStress == 1)
        cond41 = (cond4 & (tmax <= self.Tmax_lo))
        self.Kst['PolH'][cond41] = 1
        cond42 = (cond4 & (tmax >= self.Tmax_up))
        self.Kst['PolH'][cond42] = 0
        cond43 = (cond4 & np.logical_not(cond41 | cond42))
        Trel = (tmax - self.Tmax_lo) / (self.Tmax_up - self.Tmax_lo)
        self.Kst['PolH'][cond43] = ((KsPol_up * KsPol_lo) / (KsPol_lo + (KsPol_up - KsPol_lo) * np.exp(-self.fshape_b * (1 - Trel))))[cond43]

        # Calculate effects of cold stress on pollination
        cond5 = (self.PolColdStress == 0)
        self.Kst['PolC'][cond5] = 1
        cond6 = (self.PolColdStress == 1)
        cond61 = (cond6 & (tmin >= self.Tmin_up))
        self.Kst['PolC'][cond61] = 1
        cond62 = (cond6 & (tmin <= self.Tmin_lo))
        self.Kst['PolC'][cond62] = 0
        Trel = (self.Tmin_up - tmin) / (self.Tmin_up - self.Tmin_lo)
        self.Kst['PolC'][cond62] = ((KsPol_up * KsPol_lo) / (KsPol_lo + (KsPol_up - KsPol_lo) * np.exp(-self.fshape_b * (1 - Trel))))[cond62]
        
    def biomass_accumulation(self, meteo):
        """Function to calculate biomass accumulation"""

        et0 = meteo.referencePotET[None,:,:] * np.ones((nr))[:,None,None]
        self.temperature_stress(meteo)

        # Get time for harvest index build-up
        HIt = self.DAP - self.DelayedCDs - self.HIstartCD - 1

        fswitch = np.zeros((self.nRotation, self.nLat, self.nLon))
        WPadj = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond1 = (growing_season_index & (((self.CropType == 2) | (self.CropType == 3)) & (self.HIref > 0)))

        # Adjust WP for reproductive stage
        cond11 = (cond1 & (self.Determinant == 1))
        fswitch[cond11] = (self.PctLagPhase / 100)[cond11]
        cond12 = (cond1 & np.logical_not(cond11))
        cond121 = (cond12 < (self.YldFormCD / 3))
        fswitch[cond121] = (HIt / (self.YldFormCD / 3))[cond121]
        cond122 = (cond12 & np.logical_not(cond121))
        fswitch[cond122] = 1
        WPadj[cond1] = (self.WP * (1 - (1 - self.WPy / 100) * fswitch))[cond1]
        cond2 = (growing_season_index & np.logical_not(cond1))
        WPadj[cond2] = self.WP[cond2]

        # Adjust WP for CO2 effects **TODO**
        # WPadj *= self.fCO2

        # Calculate biomass accumulation on current day
        dB_NS = WPadj * (self.TrPot0 / et0) * self.Kst['Bio']  # TODO: check correct TrPot is being used
        dB = WPadj * (self.TrAct / et0) * self.Kst['Bio']  # TODO: check correct TrAct is being used

        # Update biomass accumulation
        self.B += dB
        self.B_NS += dB_NS

        # No biomass accumulation outside growing season
        self.B[np.logical_not(growing_season_index)] = 0
        self.B_NS[np.logical_not(growing_season_index)] = 0

    def harvest_index_adj_pre_anthesis(self):
        """Function to calculate adjustment to harvest index for 
        pre-anthesis water stress
        """
        # Calculate adjustment
        Br = self.B / self.B_NS
        Br_range = np.log(self.dHI_pre) / 5.62
        Br_upp = 1
        Br_low = 1 - Br_range
        Br_top = Br_upp - (Br_range / 3)

        # Get biomass ratio
        ratio_low = (Br - Br_low) / (Br_top - Br_low)
        ratio_upp = (Br - Br_top) / (Br_upp - Br_top)

        # Calculate adjustment factor
        self.Fpre = np.ones((self.nRotation, self.nLat, self.nLon))
        cond1 = (Br >= Br_low) & (Br < Br_top)
        self.Fpre[cond1] = (1 + (((1 + np.sin((1.5 - ratio_low) * np.pi)) / 2) * (self.dHI_pre / 100)))[cond1]
        cond2 = (((Br > Br_top) & (Br <= Br_upp)) & np.logical_not(cond1))
        self.Fpre[cond2] = (1 + (((1 + np.sin((0.5 + ratio_upp) * np.pi)) / 2) * (self.dHI_pre / 100)))[cond2]

        # No green canopy left at start of flowering so no harvestable crop
        # will develop
        self.Fpre[self.CC <= 0.01] = 0
        
    def harvest_index_adj_pollination(self, HIt):
        """Function to calculate adjustment to harvest index for 
        failure of pollination due to water or temperature stress
        """
        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))
        FracFlow = arr_zeros
        t1 = arr_zeros
        t2 = arr_zeros
        F1 = arr_zeros
        F2 = arr_zeros
        F = arr_zeros

        # Fractional flowering on previous day
        cond1 = (HIt > 0)
        t1[cond1] = HIt[cond1] - 1
        cond11 = (t1 > 0)
        t1pct = 100 * (t1 / self.FloweringCD)
        t1pct = np.clip(t1pct, 0, 100)
        F1[cond11] = (0.00558 * np.exp(0.63 * np.log(t1pct)) - (0.000969 * t1pct) - 0.00383)[cond11]
        F1 = np.clip(F1, 0, None)

        # Fractional flowering on current day
        t2[cond1] = HIt[cond1]
        cond12 = (t2 > 0)
        t2pct = 100 * (t2 / self.FloweringCD)
        t2pct = np.clip(t2pct, 0, 100)
        F2[cond12] = (0.00558 * np.exp(0.63 * np.log(t2pct)) - (0.000969 * t2pct) - 0.00383)[cond11]
        F2 = np.clip(F2, 0, None)

        # Weight values
        cond13 = (np.abs(F1 - F2) >= 0.0000001)
        F[cond13] = (100 * ((F1 + F2) / 2) / self.FloweringCD)[cond13]
        FracFlow[cond1] = F[cond1]
        
        # Calculate pollination adjustment for current day
        dFpol = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond2 = (self.CC >= self.CCmin)
        Ks = np.minimum(self.Ksw['Pol'], self.Kst['PolC'], self.Kst['PolH'])
        dFpol[cond2] = (Ks * FracFlow * (1 + (self.exc / 100)))[cond2]

        # Calculate pollination adjustment to date
        self.Fpol += dFpol
        self.Fpol = np.clip(self.Fpol, None, 1)

    def harvest_index_adj_post_anthesis(self):
        """Function to calculate adjustment to harvest index for 
        post-anthesis water stress
        """
        # 1 Adjustment for leaf expansion
        tmax1 = self.CanopyDevEndCD - self.HIstartCD
        DAP = self.DAP - self.DelayedCDs
        cond1 = (
            (DAP <= (self.CanopyDevEndCD + 1))
            & (tmax1 > 0)
            & (self.Fpre > 0.99)
            & (self.CC > 0.001)
            & (self.a_HI > 0))
        dCor = (1 + (1 - self.Ksw['Exp']) / self.a_HI)
        self.sCor1[cond1] += (dCor / tmax1)[cond1]
        DayCor = (DAP - 1 - self.HIstartCD)
        self.fpost_upp[cond1] = ((tmax1 / DayCor) * self.sCor1)[cond1]

        # 2 Adjustment for stomatal closure
        tmax2 = self.YldFormCD
        cond2 = (
            (DAP <= (self.HIendCD + 1))
            & (tmax2 > 0)
            & (self.Fpre > 0.99)
            & (self.CC > 0.001)
            & (self.b_HI > 0))
        dCor = ((np.exp(0.1 * np.log(self.Ksw['Sto']))) * (1 - (1 - self.Ksw['Sto']) / self.b_HI))
        self.sCor2[cond2] += (dCor / tmax2)[cond2]
        DayCor = (DAP - 1 - self.HIstartCD)
        self.fpost_dwn[cond2] = ((tmax2 / DayCor) * self.sCor2)[cond2]

        # Determine total multiplier
        cond3 = (tmax1 == 0) & (tmax2 == 0)
        self.Fpost[cond3] = 1
        cond4 = np.logical_not(cond3)
        cond41 = (cond4 & (tmax2 == 0))
        self.Fpost[cond41] = self.fpost_upp[cond41]
        cond42 = (cond4 & (tmax1 <= tmax2) & np.logical_not(cond41))
        self.Fpost[cond42] = (
            self.fpost_dwn
            * (((tmax1 * self.fpost_upp) + (tmax2 - tmax1)) / tmax2))[cond42]
        cond43 = (cond4 & np.logical_not(cond41 | cond42))
        self.Fpost[cond43] = (
            self.fpost_upp
            * (((tmax2 * self.fpost_dwn) + (tmax1 - tmax2)) / tmax2))[cond43]

    def harvest_index(self, meteo):
        """Function to simulate build up of harvest index"""
        HIadj = np.zeros((self.nRotation, self.nLat, self.nLon))
        
        # Calculate stresses, after updating root zone water content
        self.root_zone_water()
        self.water_stress(meteo, beta = True)

        # Calculate temperature stress
        self.temperature_stress(meteo)

        # Get reference harvest index on current day
        HIi = self.HIref

        # Get time for harvest index build up
        HIt = self.DAP - self.DelayedCDs - self.HIstartCD - 1

        # Calculate harvest index
        cond1 = (growing_season_index & self.YieldForm & (HIt > 0))

        # Root/tuber or fruit/grain crops
        cond11 = (cond1 & ((self.CropType == 2) | (self.CropType == 3)))

        # Determine adjustment for water stress before anthesis
        cond111 = (cond11 & np.logical_not(self.PreAdj))
        self.PreAdj[cond111] = True
        self.harvest_index_adj_pre_anthesis()

        # Adjustment only for fruit/grain crops
        HImax = np.zeros((self.nRotation, self.nLat, self.nLon))  # TODO: is this in the right place?
        cond112 = (cond11 & (self.CropType == 3))
        cond1121 = (cond112 & ((HIt > 0) & (HIt <= self.FloweringCD)))
        self.harvest_index_adj_pollination(HIt)
        HImax[cond112] = (self.Fpol * self.HI0)[cond112]
        cond113 = (cond11 & np.logical_not(cond112))
        HImax[cond113] = self.HI0[cond113]

        # Determine adjustments for post-anthesis water stress
        cond114 = (cond11 & (HIt > 0))
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
        cond12 = (cond1 & np.logical_not(cond11))
        HIadj[cond12] = HIi[cond12]

        # Otherwise no build-up of harvest index if outside yield formation
        # period
        HIi = self.HI
        HIadj = self.HIadj

        # Store final values for current time step
        self.HI = HIi
        self.HIadj = HIadj

        # No harvestable crop outside of a growing season
        self.HI[np.logical_not(growing_season_index)] = 0
        self.HIadj[np.logical_not(growing_season_index)] = 0
                
