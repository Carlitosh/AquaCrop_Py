#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
import shutil
import sys
import math
import gc

import pcraster as pcr

import virtualOS as vos
import meteo
import groundwater
import landSurface

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

        # number of soil layers/compartments
        self.nLayer = int(configuration.soilOptions['nLayer'])
        self.nComp = int(configuration.soilOptions['nComp'])
        
        # prepare sub-models
        self.createSubmodels(initialState)

    @property
    def configuration(self):
        return self._configuration

    def createSubmodels(self, initialState):
        self.meteo = meteo.Meteo(self._configuration,self.landmask,initialState)
        self.groundwater = groundwater.Groundwater(self._configuration,self.landmask,initialState)
        self.landSurface = landSurface.LandSurface(self._configuration,self._modelTime,self.meteo,self.groundwater,self.landmask,initialState)
        # TODO: identify other submodels

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

    def read_forcings(self):
        logger.info("Reading forcings for time %s", self._modelTime)
        self.meteo.read_forcings(self._modelTime)
        self.groundwater.read_forcings(self._modelTime)

    def update(self):
        logger.info("Updating model for time %s", self._modelTime)
        self.meteo.update(self._modelTime)
        self.landSurface.update(self.meteo, self.groundwater, self._modelTime)
        

    
    
