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
import SoilParameters as soilParams

class RootZoneWater(object):

    def __init__(self, iniItems, landmask):
        
        self.cloneMap = iniItems.cloneMap
        self.landmask = landmask
        
        attr = vos.getMapAttributesALL(self.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])

        # Get soil and rotation characteristics from netCDF
        soil = soilParams.SoilAndTopoParameters(iniItems, self.landmask)
        soil.read()

        rotationParameterFileNC = iniItems.rotationOptions['rotationParameterNC']
        CropSequence = vos.netcdf2PCRobjCloneWithoutTime(
            rotationParameterFileNC,
            'CropSequence',
            cloneMapFileName=self.cloneMap)

        self.nComp = soil.nComp        
        self.nRotation = CropSequence.shape[1]
        
    def update(self, th, Zroot):
        """Function to update actual and total available water in the 
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
        rootdepth = np.maximum(Zmin, self.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        comp_sto = ((dzsum - dz) < rootdepth)
        
        # Fraction of compartment covered by root zone (zero in compartments
        # NOT covered by the root zone)
        factor = np.clip(rootdepth / dzsum, 0, 1)
        factor[np.logical_not(comp_sto)] = 0

        # Actual water storage in root zone (mm)
        Wr_comp = factor * 1000 * self.th * dz
        Wr = np.sum(Wr_comp, axis=0)
        Wr[Wr < 0] = 0

        # Water storage in root zone at saturation (mm)
        WrS_comp = factor * 1000 * th_s * dz
        WrS = np.sum(WrS_comp, axis=0)
        
        # Water storage in root zone at field capacity (mm)
        WrFC_comp = factor * 1000 * th_fc * dz
        WrFC = np.sum(WrFC_comp, axis=0)

        # Water storage in root zone at permanent wilting point (mm)
        WrWP_comp = factor * 1000 * th_wp * dz
        WrWP = np.sum(WrWP_comp, axis=0)
        
        # Water storage in root zone at air dry (mm)
        WrDry_comp = factor * 1000 * th_dry * dz
        WrDry = np.sum(WrDry_comp, axis=0)

        # Water storage in root zone at aeration stress threshold (mm)
        WrAer_comp = factor * 1000 * (th_s - (self.crop.Aer / 100)) * dz
        WrAer = np.sum(WrAer_comp, axis=0)

        # Convert depths to m3/m3
        self.thRZ_Act = Wr / (rootdepth * 1000)
        self.thRZ_Sat = WrS / (rootdepth * 1000)
        self.thRZ_Fc  = WrFC / (rootdepth * 1000)
        self.thRZ_Wp  = WrWP / (rootdepth * 1000)
        self.thRZ_Dry = WrDry / (rootdepth * 1000)
        self.thRZ_Aer = WrAer / (rootdepth * 1000)

        # Calculate total available water and root zone depletion
        self.TAW = np.clip((WrFC - WrWP), 0, None)
        self.Dr = np.clip((WrFC - Wr), 0, None)
        
