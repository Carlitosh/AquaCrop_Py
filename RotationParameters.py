#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# AquaCrop
#
import os
import numpy as np
import pcraster as pcr
import virtualOS as vos
import netCDF4 as nc
import datetime as datetime
import calendar as calendar

import CropParameters as cropParams

class Rotation(object):

    def __init__(self, iniItems, landmask):
        
        self.cloneMap = iniItems.cloneMap
        self.tmpDir   = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']
        self.landmask = landmask

        attr = vos.getMapAttributesALL(self.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])

        # Crop sequence
        rotationParameterFileNC = iniItems.rotationOptions['rotationParameterNC']
        CropSequence = vos.netcdf2PCRobjCloneWithoutTime(
            rotationParameterFileNC,
            'CropSequence',
            cloneMapFileName=self.cloneMap)

        CropSequence = (
            CropSequence[:,:,None,None]
            * np.ones((self.nLat, self.nLon))[None,None,:,:])
        self.CropSequence = CropSequence.astype(bool)
        
        self.nCrop = rotation.CropSequence.shape[0]
        self.nRotation = self.CropSequence.shape[1]

        # Parameters (possibly do this in AquaCrop, and pass as arguments to __init__??)
        # self.soil = soil_parameters
        # self.crop = crop_parameters
        # self.irrigation_mgmt = irrigation_mgmt_parameters
        # self.field_mgmt = field_mgmt_parameters
        
        # Soil parameters
        self.soil = soilParams.SoilAndTopoParameters(iniItems, self.landmask)
        self.soil.read()

        # Crop parameters
        self.crop = cropParams.CropParameters(iniItems, landmask)
        self.crop.read_crop_parameters()
        self.crop.compute_variables()
        arr_ones = np.ones((self.nRotation))[None,:,None,None]
        self.PlantingDate = self.crop.PlantingDate[:,None,:,:] * arr_ones
        self.HarvestDate = self.crop.HarvestDate[:,None,:,:] * arr_ones

        # Irrigation management parameters
        self.irrigation_mgmt = irriMgmtParams.IrrigationMgmtParameters(iniItems, landmask)
        self.irrigation_mgmt.read_irrigation_mgmt_parameters()

        # Field management parameters
        self.field_mgmt = fieldMgmtParams.FieldMgmtParameters(iniItems, landmask)
        self.field_mgmt.read_field_mgmt_parameters()
        
        # Initialize season counter
        self.SeasonCounter = np.zeros((self.nCrop, self.nRotation, self.nLat, self.nLon))

        # Initialize days after planting counter
        self.DAP = np.zeros((self.nRotation, self.nLat, self.nLon))

        # Initialize cumulative GDD
        self.GDDcum = np.zeros((self.nRotation, self.nLat, self.nLon))
        
        # Initialize growing season
        self.GrowingSeason = np.full((self.nCrop, self.nRotation, self.nLat, self.nLon), False)
        
        # Initialize variables
        for nm in self.crop.parameter_names:
            vars(self)[nm] = np.zeros((self.nRotation, self.nLat, self.nLon))

        # Initialize...
        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))

        self.PreIrr = arr_zeros
        self.thRZ = {'Act': arr_zeros,
                     'Sat': arr_zeros,
                     'Fc': arr_zeros,
                     'Wp': arr_zeros,
                     'Dry': arr_zeros,
                     'Aer': arr_zeros}
        self.TAW = arr_zeros
        self.Dr = arr_zeros
        self.GrowthStage = arr_zeros
        self.Irr = arr_zeros
        self.IrrCum = arr_zeros
        # ...
        self.DelayedGDDs = arr_zeros
        self.DelayedCDs = arr_zeros
        self.GrowthStage = arr_zeros
        self.Ksw = {'Exp': arr_zeros,
                    'Sto': arr_zeros,
                    'Sen': arr_zeros,
                    'Pol': arr_zeros,
                    'StoLin': arr_zeros}
        
        self.Kst = {'Bio': arr_zeros,
                    'PolH': arr_zeros,
                    'PolC': arr_zeros}
        
    def growing_degree_day(self, meteo):
        """Function to calculate number of growing degree days on 
        current day
        """
        tmax = meteo.tmax[None,:,:] * np.ones((self.nCrop))[:,None,None]
        tmin = meteo.tmin[None,:,:] * np.ones((self.nCrop))[:,None,None]
        
        if self.GDDmethod == 1:
            tmean = ((tmax + tmin) / 2)
            tmean = np.clip(tmean, self.Tbase, self.Tupp)
        elif self.GDDmethod == 2:
            tmax = np.clip(tmax, self.Tbase, self.Tupp)
            tmin = np.clip(tmin, self.Tbase, self.Tupp)
            tmean = ((tmax + tmin) / 2)
        elif self.GDDmethod == 3:
            tmax = np.clip(tmax, self.Tbase, self.Tupp)
            tmin = np.clip(tmin, None, self.Tupp)
            tmean = np.clip(tmean, self.Tbase, None)
            
        self.GDD = (tmean - self.Tbase)

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
        PreIrr_comp = ((thCrit - self.th) * 1000 * dz)
        PreIrr_comp[np.logical_not(cond1)] = 0
        self.PreIrr_comp = PreIrr_comp
        self.PreIrr = np.sum(PreIrr_comp, axis=0)
        self.th[cond1] = thCrit[cond1]
        
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
        if self.crop.CalendarType == 1:
            tAdj = (self.DAP - self.DelayedCDs)
        elif self.crop.CalendarType == 2:
            tAdj = (self.GDDcum - self.DelayedGDDs)
            
        # Calculate root expansion
        Zini = self.Zmin * (self.PctZmin / 100)
        t0 = np.round(self.Emergence / 2)
        tmax = self.MaxRooting
        if self.crop.CalendarType == 1:
            tOld = (tAdj - 1)
        elif self.crop.CalendarType == 2:
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
        TrAct0 = np.zeros((self.nRotation, self.nLat, self.nLon))
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
        self.SurfaceStorage[cond101] -= (fSub * self.TrPot0)[cond101]
        TrAct0[cond101] = (fSub * self.TrPot0)[cond101]

        # Otherwise there is no transpiration from surface storage
        cond102 = (cond10 & np.logical_not(cond101))
        TrAct0[cond102] = 0

        # More water can be extracted from soil profile for transpiration
        cond103 = (cond10 & (TrAct0 < (fSub * self.TrPot0)))
        TrPot[cond103] = ((fSub * self.TrPot0) - TrAct0)[cond103]

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
        self.TrAct += TrAct0

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
        
    def update(self, meteo, currTimeStep):
        """Function to update parameters for current crop grown as well 
        as counters pertaining to crop growth
        """
        # Update season counter for crops planted on current day
        cond1 = (self.CropSequence & (currTimeStep.doy == self.PlantingDate))
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
            vars(self)[var] = param[crop_index,I,J,K]

        for nm in self.irrigation_mgmt.parameter_names:
            vars(self)[var] = getattr(self.irrigation_mgmt, nm)[crop_index,I,J,K]

        for nm in self.field_mgmt.parameter_names:
            vars(self)[var] = getattr(self.field_mgmt, nm)[crop_index,I,J,K]
