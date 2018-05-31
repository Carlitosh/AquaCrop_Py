#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
from pcraster.framework import *
import pcraster as pcr
import numpy as np

import logging
logger = logging.getLogger(__name__)

import VirtualOS as vos
# from ncConverter import *
# import ETPFunctions as refPotET

import logging
logger = logging.getLogger(__name__)

class CO2(object):

    def __init__(self, iniItems, landmask, initialState = None):
        object.__init__(self)

        self.cloneMap = iniItems.cloneMap
        self.tmpDir = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']
        
        # landmask/area of interest
        self.landmask = landmask
        if iniItems.globalOptions['landmask'] != "None":
           self.landmask = vos.readPCRmapClone(iniItems.globalOptions['landmask'],\
                                               self.cloneMap,self.tmpDir,self.inputDir)

        # reference concentration
        self.RefConc = 369.41
        
        self.co2FileNC = iniItems.carbonDioxideOptions['carbonDioxideNC']

        # variable names      
        self.co2VarName = 'co2' 
        if 'co2VariableName' in iniItems.carbonDioxideOptions:
            self.co2VarName = carbonDioxideOptions['co2VariableName']

        self.co2_set_per_year  = False
        
        # make the iniItems available for the other modules:
        self.iniItems = iniItems

    def update(self,currTimeStep):
        #TODO: calculate  referencePotET
        pass

    def read_forcings(self, currTimeStep):

        # TODO: only read data if the year has changed
        
        date = '%04i-%02i-%02i' %(currTimeStep.year, 1, 1)
        self.conc = vos.netcdf2PCRobjClone(self.co2FileNC,
                                           self.co2VarName,
                                           date,
                                           useDoy = None,
                                           cloneMapFileName = self.cloneMap,
                                           LatitudeLongitude = True)
