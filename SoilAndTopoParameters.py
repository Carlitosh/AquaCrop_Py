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

class SoilAndTopoParameters(object):

    def __init__(self, iniItems, landmask):
        object.__init__(self)

        # cloneMap, tmpDir, inputDir based on the configuration/setting given in the ini/configuration file
        self.cloneMap = iniItems.cloneMap
        self.tmpDir   = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']
        self.landmask = landmask  # TODO: do something with this - mask out cells outside the model area?

        attr = vos.getMapAttributesALL(self.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])
        
        # characteristics of soil layers/compartments
        self.nLayer = int(iniItems.soilOptions['nLayer'])
        self.zLayer = iniItems.soilOptions['zLayer'].split(',')
        self.zLayer = np.array(map(np.float32, self.zLayer))

        self.nComp = int(iniItems.soilOptions['nComp'])
        self.dz = iniItems.soilOptions['dz'].split(',')
        self.dz = np.array(map(np.float32, self.dz))
        self.dzsum = np.round(100 * np.cumsum(self.dz)) / 100

        # soil parameter input file
        self.soilAndTopoFileNC = iniItems.soilOptions['soilAndTopoNC']

    def read(self):
		
        self.readTopo()
        self.readSoil()

    def readTopo(self):
        pass

    def readSoil(self):

        soilParams = ['ksat','th_s','th_fc','th_wp','CalcSHP','EvapZsurf','EvapZmin','EvapZmax','Kex','fevap','fWrelExp','fwcc','AdjREW','REW','AdjCN','CN','zCN','zGerm','zRes','fshape_cr']
        for var in soilParams:
            vars(self)[var] = vos.netcdf2PCRobjCloneWithoutTime(self.soilAndTopoFileNC,\
                                                                var,\
                                                                cloneMapFileName=self.cloneMap)

        self.nRotation = self.ksat.shape[1]
        
        # map layers to compartments - the result is a 1D array with length
        # equal to nComp where the value of each element is the index of the
        # corresponding layer. For the time being we use the layer in which
        # the midpoint of each compartment is located.
        #
        # some data to try the command yourself:
        # 
        # from numpy import np
        # zMid = np.array((0.05,0.15,0.25,0.375,0.525,0.7,0.9,1.125,1.375,1.625,1.875,2.15))
        # zLayerTop = np.array((0,0.5))
        # nLayer = 2
        zBot = np.cumsum(self.dz)
        zTop = zBot - self.dz
        zMid = (zTop + zBot) / 2
        zLayerBot = np.cumsum(self.zLayer)
        zLayerTop = zLayerBot - self.zLayer        
        self.layerIndex = np.sum(((zMid[:,None] * np.ones((self.nLayer))[None,:]) > zLayerTop), axis=1) - 1

        # # transform certain soil properties to (ncomp, nrotation, nlat, nlon)
        # soilParams = ['ksat','th_s','th_fc','th_wp']
        # for var in soilParams:
        #     vars(self)[var] = vars(self)[var][self.layerIndex,:,:,:]

        # The following is adapted from AOS_ComputeVariables.m, lines 129-139
        # "Calculate drainage characteristic (tau)
        self.tau = 0.0866 * (self.ksat ** 0.35)
        self.tau = np.round(self.tau * 100) / 100
        self.tau[self.tau > 1] = 1
        self.tau[self.tau < 0] = 0

        # The following is adapted from AOS_ComputeVariables.m, lines 25
        self.th_dry = self.th_wp / 2

        # The following is adapted from AOS_ComputeVariables.m, lines 147-151
        # "Calculate upper and lower curve numbers"
        self.CNbot = np.round(1.4 * (np.exp(-14 * np.log(10))) + (0.507 * self.CN) - (0.00374 * self.CN ** 2) + (0.0000867 * self.CN ** 3))
        self.CNtop = np.round(5.6 * (np.exp(-14 * np.log(10))) + (2.33 * self.CN) - (0.0209 * self.CN ** 2) + (0.000076 * self.CN ** 3))

    def compute_capillary_rise_parameters(self):
        # Function adapted from AOS_ComputeVariables.m, lines 60-127

        self.aCR = np.zeros((self.nLayer, self.nRotation, self.nLat, self.nLon))
        self.bCR = np.zeros((self.nLayer, self.nRotation, self.nLat, self.nLon))

        # "Sandy soil class"
        cond1 = (self.th_wp >= 0.04) & (self.th_wp <= 0.15) & (self.th_fc >= 0.09) & (self.th_fc <= 0.28) & (self.th_s >= 0.32) & (self.th_s <= 0.51)
        cond11 = (cond1 & (self.ksat >= 200) & (self.ksat <= 2000))
        cond12 = (cond1 & (self.ksat < 200))
        cond13 = (cond1 & (self.ksat > 2000))
        self.aCR[cond11] = (-0.3112 - (self.ksat * (10 ** -5)))[cond11]
        self.bCR[cond11] = (-1.4936 + (0.2416 * np.log(self.ksat)))[cond11]
        self.aCR[cond12] = (-0.3112 - (200 * (10 ** -5)))
        self.bCR[cond12] = (-1.4936 + (0.2416 * np.log(200)))
        self.aCR[cond13] = (-0.3112 - (2000 * (10 ** -5)))
        self.bCR[cond13] = (-1.4936 + (0.2416 * np.log(2000)))

        # "Loamy soil class"
        cond2 = (self.th_wp >= 0.06) & (self.th_wp <= 0.20) & (self.th_fc >= 0.23) & (self.th_fc <= 0.42) & (self.th_s >= 0.42) & (self.th_s <= 0.55)
        cond21 = (cond2 & (self.ksat >= 100) & (self.ksat <= 750))
        cond22 = (cond2 & (self.ksat < 100))
        cond23 = (cond2 & (self.ksat > 750))
        self.aCR[cond21] = (-0.4986 + (9 * (10 ** -5) * self.ksat))[cond21]
        self.bCR[cond21] = (-2.1320 + (0.4778 * np.log(self.ksat)))[cond21]
        self.aCR[cond22] = (-0.4986 + (9 * (10 ** -5) * 100))
        self.bCR[cond22] = (-2.1320 + (0.4778 * np.log(100)))
        self.aCR[cond23] = (-0.4986 + (9 * (10 ** -5) * 750))
        self.bCR[cond23] = (-2.1320 + (0.4778 * np.log(750)))

        # "Sandy clayey soil class"
        cond3 = (self.th_wp >= 0.16) & (self.th_wp <= 0.34) & (self.th_fc >= 0.25) & (self.th_fc <= 0.45) & (self.th_s >= 0.40) & (self.th_s <= 0.53)
        cond31 = (cond3 & (self.ksat >= 5) & (self.ksat <= 150))
        cond32 = (cond3 & (self.ksat < 5))
        cond33 = (cond3 & (self.ksat > 150))
        self.aCR[cond31] = (-0.5677 - (4 * (10 ** -5) * self.ksat))[cond31]
        self.bCR[cond31] = (-3.7189 + (0.5922 * np.log(self.ksat)))[cond31]
        self.aCR[cond32] = (-0.5677 - (4 * (10 ** -5) * 5))
        self.bCR[cond32] = (-3.7189 + (0.5922 * np.log(5)))
        self.aCR[cond33] = (-0.5677 - (4 * (10 ** -5) * 150))
        self.bCR[cond33] = (-3.7189 + (0.5922 * np.log(150)))

        # "Silty clayey soil class"
        cond4 = (self.th_wp >= 0.20) & (self.th_wp <= 0.42) & (self.th_fc >= 0.40) & (self.th_fc <= 0.58) & (self.th_s >= 0.49) & (self.th_s <= 0.58)
        cond41 = (cond4 & (self.ksat >= 1) & (self.ksat <= 150))
        cond42 = (cond4 & (self.ksat < 1))
        cond43 = (cond4 & (self.ksat > 150))
        self.aCR[cond41] = (-0.6366 + (8 * (10 ** -4) * self.ksat))[cond41]
        self.bCR[cond41] = (-1.9165 + (0.7063 * np.log(self.ksat)))[cond41]
        self.aCR[cond42] = (-0.6366 + (8 * (10 ** -4) * 1))
        self.bCR[cond42] = (-1.9165 + (0.7063 * np.log(1)))
        self.aCR[cond43] = (-0.6366 + (8 * (10 ** -4) * 150))
        self.bCR[cond43] = (-1.9165 + (0.7063 * np.log(150)))

