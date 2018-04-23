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

class SoilWaterBalance(object):

    def __init__(self, configuration, currTimeStep, initialState = None):

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

        # Soil parameters
        self.soil = soilParams.SoilAndTopoParameters(iniItems, self.landmask)
        self.soil.read()

        # Compute capillary rise parameters if water table is modelled
        if groundwater.WaterTable:
            self.soil.compute_capillary_rise_parameters()

        # If groundwater is present, this function calculates adjusted FC
        self.soil.check_groundwater_table(groundwater)
        
        # Define variables
        # TODO: take out of function
        self.set_initial_conditions()

        # list of state variables
        self.state_vars = ['AgeDays','AgeDays_NS','AerDays','IrrCum','IrrNetCum',
                           'DaySubmerged','Epot','Tpot','WTinSoil','AerDaysComp',
                           'TrRatio','SurfaceStorage','SurfaceStorageIni','th',
                           'Wsurf','EvapZ','Stage2','Wstage2']
        
    def set_initial_conditions(self):

        # TODO: set initial conditions
        # self.init_conditions = initCond.InitialConditions(self._configuration, self._modelTime, self.crop)

        # Groundwater
        self.GwIn = arr_zeros
        self.zGW = arr_zeros

        # Drainage
        self.FluxOut = arr_zeros
        self.DeepPerc = arr_zeros

        # Infiltration
        self.Runoff = arr_zeros
        self.Infl = arr_zeros

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
        self.Irr = arr_zeros
        self.PreIrr = arr_zeros
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

        # Aeration stress
        self.AerDaysComp = np.zeros((self.nComp, self.nRotation, self.nLat, self.nLon))
                
        # Transpiration
        self.Ksa_Aer = arr_zeros
        self.TrPot0 = arr_zeros
        self.TrPot_NS = arr_zeros
        self.TrAct = arr_zeros
        self.TrAct0 = arr_zeros
        self.Tpot = arr_zeros

    def getState(self):
        result = {}
        for var in self.state_vars:
            result[var] = vars(self)[var]
        return result

    def getPseudoState(self):
        pass

    def update(self, meteo, currTimeStep):
        pass
    
    def pre_irrigation(self, crop):
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
        rootdepth = np.maximum(crop.Zmin, crop.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        thCrit = (th_wp + ((crop.NetIrrSMT / 100) * (th_fc - th_wp)))

        # Conditions for applying pre-irrigation
        cond1 = ((crop.IrrMethod == 4) & (crop.DAP == 1) & ((dzsum - dz) < rootdepth) & (self.th < thCrit))

        # Update pre-irrigation and root zone water content (mm)
        PreIrr_req = ((thCrit - self.th) * 1000 * dz)
        PreIrr_req[np.logical_not(cond1)] = 0
        # self.PreIrr_req = PreIrr_req
        self.PreIrr = np.sum(PreIrr_req, axis=0)

    def drainage(self):
        """Function to redistribute stored soil water"""
        
        # Expand soil properties to compartments
        ksat = self.soil.ksat[self.soil.layerIndex,:]
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_s = self.soil.th_s[self.soil.layerIndex,:]
        tau = self.soil.tau[self.soil.layerIndex,:]
        
        # Preallocate arrays        
        thnew = self.th
        drainsum = np.zeros((self.nRotation, self.nLat, self.nLon))

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
            cond643 = (
                cond64
                & np.logical_not(cond641 | cond642))

            # Drainage and cumulative drainage are zero as water content has not
            # risen above field capacity in the present compartment
            drainsum[cond643] = 0
                
            # ##################################################################
            # Increase water content in compartment ii with cumulative drainage
            # from above
            cond65 = (np.logical_not(drainability) & (thX > th_s[comp,:]))
            thnew[comp,:][cond65] = (
                self.th[comp,:]
                + (drainsum / (1000 * self.soil.dz[comp])))[cond65]

            # Check new water content against hydraulic properties of soil layer
            cond651 = (cond65 & (thnew[comp,:] <= th_s[comp,:]))

            # Calculate new drainage ability
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

            cond6512 = (cond651 & (np.logical_not(cond6511)))
            drainsum[cond6512] = 0
            
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

    def rainfall_partition(self, crop, meteo):
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
        cond1 = ((crop.Bunds == 0) | (crop.zBund < 0.001))
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
        CN[cond11] = np.round(self.soil.CNbot + (self.soil.CNtop - self.soil.CNbot) * wet_top)[cond11]

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

    def root_zone_water(self, soil):
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
        Wr_comp = factor * 1000 * self.th * dz
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

    def irrigation(self, crop, meteo):
        """Function to get irrigation depth for the current day"""

        # Expand soil properties to compartments
        th_s   = self.soil.th_s[self.soil.layerIndex,:]
        th_fc  = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp  = self.soil.th_wp[self.soil.layerIndex,:]
        th_dry = self.soil.th_dry[self.soil.layerIndex,:]

        SMT = np.concatenate((crop.SMT1[None,:],
                              crop.SMT2[None,:],
                              crop.SMT3[None,:],
                              crop.SMT4[None,:]), axis=0)

        # Calculate root zone water content and depletion
        self.root_zone_water()

        # Determine adjustment for inflows and outflows on current day
        cond1 = (self.thRZ['Act'] > self.thRZ['Fc'])
        rootdepth = np.zeros((self.nRotation, self.nLat, self.nLon))
        rootdepth[cond1] = np.maximum(crop.Zmin, crop.Zroot)[cond1]
        AbvFc = np.zeros((self.nRotation, self.nLat, self.nLon))
        AbvFc[cond1] = ((self.thRZ['Act'] - self.thRZ['Fc']) * 1000 * rootdepth)[cond1]
        WCadj = self.Tpot + self.Epot - meteo.precipitation + self.Runoff - AbvFc

        # Determine irrigation depth (mm/day) to be applied

        # TODO: this should be placed elsewhere, GrowthStage should be a crop parameter
        # Update growth stage if it is first day of a growing season
        cond2 = (growing_season_index & (crop.DAP == 1))
        self.GrowthStage[cond2] = 1

        # Run irrigation depth calculation

        # No irrigation if rainfed (i.e. IrrMethod == 0)
        condX = (growing_season_index & (crop.IrrMethod == 0))
        self.Irr[condX] = 0
        
        # If irrigation is based on soil moisture, get the soil moisture
        # target for the current growth stage and determine threshold to
        # initiate irrigation
        cond3 = (growing_season_index & (crop.IrrMethod == 1))
        I,J,K = np.ogrid[:self.nRotation,:self.nLat,:self.nLon]
        growth_stage_index = crop.GrowthStage - 1
        growth_stage_index = growth_stage_index.astype(int)
        SMT = SMT[growth_stage_index,I,J,K]
        IrrThr = (1 - SMT / 100) * self.TAW
                
        # If irrigation is based on a fixed interval, get number of days in
        # growing season so far (subtract 1 so that we always irrigate first
        # on day 1 of each growing season)
        cond4 = (growing_season_index & (self.IrrMethod == 2))
        nDays = crop.DAP - 1

        # combine irrigation methods 1 and 2
        cond5 = (cond3 | cond4)
        
        # Working on a copy, adjust depletion for inflows and outflows - same
        # for both soil moisture and interval based irrigation
        Dr = self.Dr
        Dr[cond5] += WCadj[cond5]
        Dr[cond5] = np.clip(Dr, 0, None)[cond5]

        # check if conditions for irrigation method 1 or 2 are met
        cond6 = ((cond3 & (Dr > IrrThr)) | (cond4 & ((nDays % IrrInterval) == 0)))
        
        IrrReq = Dr
        EffAdj = ((100 - self.AppEff) + 100) / 100
        IrrReq *= EffAdj
        self.Irr[cond6] = np.clip(IrrReq, 0, crop.MaxIrr)[cond6]

        cond7 = (cond5 & np.logical_not(cond6))
        self.Irr[cond7] = 0

        # If irrigation is based on a pre-defined schedule then the irrigation
        # requirement for each rotation is read from a netCDF file. Note that if
        # the option 'irrScheduleFileNC' is None, then nothing will be imported
        # and the irrigation requirement will be zero
        cond8 = (growing_season_index & (crop.IrrMethod == 3))
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
        cond9 = (growing_season_index & (crop.IrrMethod == 4))
        self.Irr[cond9] = 0

        # No irrigation outside the growing season
        self.Irr[np.logical_not(growing_season_index)] = 0

        # Update cumulative irrigation counter for growing season
        self.IrrCum += self.Irr
        
    def infiltration(self, crop):
        """Function to infiltrate incoming water (rainfall and 
        irrigation)
        
        Updates Infl, DeepPerc and Runoff.
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
        self.Infl += self.Irr * (crop.AppEff / 100)  # TODO: has efficiency not already been accounted for?

        # Determine surface storage if bunds are present
        cond1 = (crop.Bunds == 1)
        cond11 = (cond1 & (crop.zBund > 0.001))
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

        # If total infiltration is zero then there is no storage or runoff
        cond112 = (cond11 & np.logical_not(cond111))
        ToStore[cond112] = 0
        RunoffIni[cond112] = 0
        
        # If bunds are not on field then infiltration is limited by saturated
        # hydraulic conductivity of top soil layer
        cond2 = (crop.Bunds == 0)
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
                 & (crop.Bunds == 1)
                 & crop.zBund > 0.001)
        self.SurfaceStorage += (Runoff - RunoffIni)

        # Limit surface storage to bund height: additional water above top of
        # bunds becomes runoff, and surface storage equals bund height
        cond51 = (cond5 & (self.SurfaceStorage > (crop.zBund * 1000)))
        Runoff[cond51] = (
            RunoffIni
            + (self.SurfaceStorage - (crop.zBund * 1000)))[cond51]
        self.SurfaceStorage[cond51] = (crop.zBund * 1000)[cond51]
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
        cond311 = (cond31 & (self.CC > self.CCxAct))
        mult[cond311] = 0
        cond312 = (cond31 & np.logical_not(cond311))
        mult[cond312] = ((self.CCxAct - self.CC) / (self.CCxAct / 2))[cond312]
        cond32 = (cond3 & np.logical_not(cond31))
        mult[cond32] = 1
        EsPot[cond3] = (EsPot * (1 - self.CCxAct * (self.soil.fwcc / 100) * mult))[cond3]
        CCxActAdj = ((1.72 * self.CCxAct) + (self.CCxAct ** 2) - 0.3 * (self.CCxAct ** 3))
        EsPotMin = np.zeros((nr, nlat, nlon))
        EsPotMin[cond3] = (self.soil.Kex * (1 - CCxActAdj) * et0)[cond3]
        EsPotMin = np.clip(EsPotMin, 0, None)

        # Line 85-89 of AOS_SoilEvaporation.m
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
        cond52 = (cond5 & (self.Mulches == 1))
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
                comp_sto = ((dzsum - dz) < self.soil.EvapZ)
                factor = 1 - ((dzsum - self.soil.EvapZ) / dz)

                # multiply by comp_sto to ensure factor is zero in compartments entirely below EvapZ
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
        DAPadj = (crop.DAP - crop.DelayedCDs) * growing_season_index
        cond5 = (growing_season_index & (DAPadj > crop.MaxCanopyCD))
        self.AgeDays[cond5] = (DAPadj - crop.MaxCanopyCD)[cond5]

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

        cond11 = (growing_season_index & np.logical_not(cond10))
        self.TrPot[cond11] = 0
        self.TrAct0[cond11] = 0
        
        # ######################################################################
        # Update potential root zone transpiration for water stress

        # Determine root zone water content
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
                # as water does not exceed wilting point
                KsComp[(cond14 & np.logical_not(cond142 | cond143))] = 0

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
