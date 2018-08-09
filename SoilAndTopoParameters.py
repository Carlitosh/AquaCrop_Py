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

    def __init__(self, SoilAndTopoParameters_variable):
        self.var = SoilAndTopoParameters_variable

    def initial(self):

        # self.var.nLayer = int(self.var._configuration.soilOptions['nLayer'])
        # self.var.zLayer = self.var._configuration.soilOptions['zLayer'].split(',')
        # self.var.zLayer = np.array(map(np.float32, self.var.zLayer))

        # self.var.nComp = int(self.var._configuration.soilOptions['nComp'])
        # self.var.dz = self.var._configuration.soilOptions['dz'].split(',')
        # self.var.dz = np.array(map(np.float64, self.var.dz))
        # self.var.dzsum = np.cumsum(self.var.dz)

        # read parameters
        self.var.soilAndTopoFileNC = self.var._configuration.soilOptions['soilAndTopoNC']
        self.read()
        
        # for convenience
        self.var.dz_xy = self.var.dz[None,:,None,None] * np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))[:,None,:,:]
        self.var.dzsum_xy = self.var.dzsum[None,:,None,None] * np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))[:,None,:,:]
        
        # # read parameters
        # self.var.soilAndTopoFileNC = self.var._configuration.soilOptions['soilAndTopoNC']
        # self.read()
        self.compute_capillary_rise_parameters()

    def read(self):		
        self.readTopo()
        self.readSoil()

    def readTopo(self):
        pass

    def readSoil(self):

        # Get layer and compartment thickness (zLayer and dz) from dimensions
        # Layers:
        zlayermid = vos.netcdfDim2NumPy(self.var.soilAndTopoFileNC,'layer')
        nlayer = zlayermid.size
        zlayer = np.zeros((nlayer))
        runtot = 0
        for layer in range(nlayer):
            zlayer[layer] = (zlayermid[layer] - runtot) * 2
            runtot = np.sum(zlayer[:(layer+1)])
        self.var.nLayer = nlayer
        self.var.zLayer = zlayer
        
        # Compartments:
        zmid = vos.netcdfDim2NumPy(self.var.soilAndTopoFileNC,'compartment')
        ncomp = zmid.size
        dz = np.zeros((ncomp))
        runtot = 0
        for comp in range(ncomp):
            dz[comp] = (zmid[comp] - runtot) * 2
            runtot = np.sum(dz[:(comp+1)])
        self.var.nComp = ncomp
        self.var.dz = dz
        self.var.dzsum = np.cumsum(dz)
        
        # These parameters have dimensions depth,lat,lon
        soilParams1 = ['ksat', 'th_s', 'th_fc', 'th_wp']
        for var in soilParams1:
            d = vos.netcdf2PCRobjCloneWithoutTime(self.var.soilAndTopoFileNC,
                                                  var,
                                                  cloneMapFileName=self.var.cloneMap)
            vars(self.var)[var] = np.broadcast_to(d, (self.var.nCrop, self.var.nLayer, self.var.nLat, self.var.nLon))

        # These parameters have dimensions lat,lon
        soilParams2 = ['CalcSHP', 'EvapZsurf','EvapZmin', 'EvapZmax', 'Kex',
                       'fevap', 'fWrelExp', 'fwcc','AdjREW', 'REW', 'AdjCN',
                       'CN', 'zCN', 'zGerm', 'zRes', 'fshape_cr']
        for var in soilParams2:
            d = vos.netcdf2PCRobjCloneWithoutTime(self.var.soilAndTopoFileNC,
                                                  var,
                                                  cloneMapFileName=self.var.cloneMap)
            d = np.broadcast_to(d, (self.var.nCrop, self.var.nLat, self.var.nLon))
            vars(self.var)[var] = np.copy(d)#np.broadcast_to(d, (self.var.nCrop, self.var.nLat, self.var.nLon))

        # map layers to compartments - the result is a 1D array with length
        # equal to nComp where the value of each element is the index of the
        # corresponding layer. For the time being we use the layer in which
        # the midpoint of each compartment is located.
        #
        # some data to try the command yourself.var:
        # 
        # from numpy import np
        # zMid = np.array((0.05,0.15,0.25,0.375,0.525,0.7,0.9,1.125,1.375,1.625,1.875,2.15))
        # zLayerTop = np.array((0,0.5))
        # nLayer = 2
        zBot = np.cumsum(self.var.dz)
        zTop = zBot - self.var.dz
        zMid = (zTop + zBot) / 2
        zLayerBot = np.cumsum(self.var.zLayer)
        zLayerTop = zLayerBot - self.var.zLayer        
        self.var.layerIndex = np.sum(((zMid[:,None] * np.ones((self.var.nLayer))[None,:]) > zLayerTop), axis=1) - 1
        
        # The following is adapted from AOS_ComputeVariables.m, lines 129-139
        # "Calculate drainage characteristic (tau)
        self.var.tau = 0.0866 * (self.var.ksat ** 0.35)
        self.var.tau = np.round(self.var.tau * 100) / 100
        self.var.tau[self.var.tau > 1] = 1
        self.var.tau[self.var.tau < 0] = 0

        # The following is adapted from AOS_ComputeVariables.m, lines 25
        self.var.th_dry = self.var.th_wp / 2

        cond = (self.var.AdjREW == 0)
        self.var.REW[cond] = (np.round((1000 * (self.var.th_fc[:,0,...] - self.var.th_dry[:,0,...]) * self.var.EvapZsurf)))[cond]
        # The following is adapted from AOS_ComputeVariables.m, lines 147-151
        # "Calculate upper and lower curve numbers"
        self.var.CNbot = np.round(1.4 * (np.exp(-14 * np.log(10))) + (0.507 * self.var.CN) - (0.00374 * self.var.CN ** 2) + (0.0000867 * self.var.CN ** 3))
        self.var.CNtop = np.round(5.6 * (np.exp(-14 * np.log(10))) + (2.33 * self.var.CN) - (0.0209 * self.var.CN ** 2) + (0.000076 * self.var.CN ** 3))
        
        # transform certain soil properties to (ncrop, ncomp, nlat, nlon)
        soil_params = ['th_s','th_fc','th_wp','th_dry','ksat','tau']
        # soil_params = ['th_s','th_fc','th_wp','th_dry','ksat','tau','fshape_cr']
        for nm in soil_params:
            newnm = nm + '_comp'
            vars(self.var)[newnm] = vars(self.var)[nm][:,self.var.layerIndex,...]

    def compute_capillary_rise_parameters(self):
        # Function adapted from AOS_ComputeVariables.m, lines 60-127

        self.var.aCR = np.zeros((self.var.nCrop, self.var.nLayer, self.var.nLat, self.var.nLon))
        self.var.bCR = np.zeros((self.var.nCrop, self.var.nLayer, self.var.nLat, self.var.nLon))

        # "Sandy soil class"
        cond1 = (self.var.th_wp >= 0.04) & (self.var.th_wp <= 0.15) & (self.var.th_fc >= 0.09) & (self.var.th_fc <= 0.28) & (self.var.th_s >= 0.32) & (self.var.th_s <= 0.51)
        cond11 = (cond1 & (self.var.ksat >= 200) & (self.var.ksat <= 2000))
        cond12 = (cond1 & (self.var.ksat < 200))
        cond13 = (cond1 & (self.var.ksat > 2000))
        self.var.aCR[cond11] = (-0.3112 - (self.var.ksat * (10 ** -5)))[cond11]
        self.var.bCR[cond11] = (-1.4936 + (0.2416 * np.log(self.var.ksat)))[cond11]
        self.var.aCR[cond12] = (-0.3112 - (200 * (10 ** -5)))
        self.var.bCR[cond12] = (-1.4936 + (0.2416 * np.log(200)))
        self.var.aCR[cond13] = (-0.3112 - (2000 * (10 ** -5)))
        self.var.bCR[cond13] = (-1.4936 + (0.2416 * np.log(2000)))

        # "Loamy soil class"
        cond2 = (self.var.th_wp >= 0.06) & (self.var.th_wp <= 0.20) & (self.var.th_fc >= 0.23) & (self.var.th_fc <= 0.42) & (self.var.th_s >= 0.42) & (self.var.th_s <= 0.55)
        cond21 = (cond2 & (self.var.ksat >= 100) & (self.var.ksat <= 750))
        cond22 = (cond2 & (self.var.ksat < 100))
        cond23 = (cond2 & (self.var.ksat > 750))
        self.var.aCR[cond21] = (-0.4986 + (9 * (10 ** -5) * self.var.ksat))[cond21]
        self.var.bCR[cond21] = (-2.1320 + (0.4778 * np.log(self.var.ksat)))[cond21]
        self.var.aCR[cond22] = (-0.4986 + (9 * (10 ** -5) * 100))
        self.var.bCR[cond22] = (-2.1320 + (0.4778 * np.log(100)))
        self.var.aCR[cond23] = (-0.4986 + (9 * (10 ** -5) * 750))
        self.var.bCR[cond23] = (-2.1320 + (0.4778 * np.log(750)))

        # "Sandy clayey soil class"
        cond3 = (self.var.th_wp >= 0.16) & (self.var.th_wp <= 0.34) & (self.var.th_fc >= 0.25) & (self.var.th_fc <= 0.45) & (self.var.th_s >= 0.40) & (self.var.th_s <= 0.53)
        cond31 = (cond3 & (self.var.ksat >= 5) & (self.var.ksat <= 150))
        cond32 = (cond3 & (self.var.ksat < 5))
        cond33 = (cond3 & (self.var.ksat > 150))
        self.var.aCR[cond31] = (-0.5677 - (4 * (10 ** -5) * self.var.ksat))[cond31]
        self.var.bCR[cond31] = (-3.7189 + (0.5922 * np.log(self.var.ksat)))[cond31]
        self.var.aCR[cond32] = (-0.5677 - (4 * (10 ** -5) * 5))
        self.var.bCR[cond32] = (-3.7189 + (0.5922 * np.log(5)))
        self.var.aCR[cond33] = (-0.5677 - (4 * (10 ** -5) * 150))
        self.var.bCR[cond33] = (-3.7189 + (0.5922 * np.log(150)))

        # "Silty clayey soil class"
        cond4 = (self.var.th_wp >= 0.20) & (self.var.th_wp <= 0.42) & (self.var.th_fc >= 0.40) & (self.var.th_fc <= 0.58) & (self.var.th_s >= 0.49) & (self.var.th_s <= 0.58)
        cond41 = (cond4 & (self.var.ksat >= 1) & (self.var.ksat <= 150))
        cond42 = (cond4 & (self.var.ksat < 1))
        cond43 = (cond4 & (self.var.ksat > 150))
        self.var.aCR[cond41] = (-0.6366 + (8 * (10 ** -4) * self.var.ksat))[cond41]
        self.var.bCR[cond41] = (-1.9165 + (0.7063 * np.log(self.var.ksat)))[cond41]
        self.var.aCR[cond42] = (-0.6366 + (8 * (10 ** -4) * 1))
        self.var.bCR[cond42] = (-1.9165 + (0.7063 * np.log(1)))
        self.var.aCR[cond43] = (-0.6366 + (8 * (10 ** -4) * 150))
        self.var.bCR[cond43] = (-1.9165 + (0.7063 * np.log(150)))

        # Expand to soil compartments
        self.var.aCR_comp = self.var.aCR[:,self.var.layerIndex,...]
        self.var.bCR_comp = self.var.bCR[:,self.var.layerIndex,...]

    def dynamic(self):
        pass
