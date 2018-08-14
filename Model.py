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

class Model(object):
    
    def __init__(self, configuration, modelTime, initialState = None):

        self._configuration = configuration
        self._modelTime = modelTime

        # clone map, land mask
        pcr.setclone(configuration.cloneMap)
        self.cloneMap = self._configuration.cloneMap
        self.landmask = vos.readPCRmapClone(configuration.globalOptions['landmask'],
                                            configuration.cloneMap,
                                            configuration.tmpDir,
                                            configuration.globalOptions['inputDir'],
                                            True)
        self.landmask = self.landmask > 0  # boolean
        
        attr = vos.getMapAttributesALL(self.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])
        
    @property
    def configuration(self):
        return self._configuration

    # def dumpState(self, outputDirectory, specific_date_string = None):

    #     if specific_date_string is None:
    #         specific_date_string = str(self._modelTime.fulldate)

    #     state = self.getState()
    #     groundwater_state = state['groundwater']
    #     for var in groundwater_state.iterItems():
    #         fn = self.outNCDir+"/"+str(var)+specific_date_string+".npy"
    #         np.save(fn, groundwater_state[var])
        
    #     soil_water_state = state['soil_water']
    #     for var in soil_water_state.iterItems():
    #         fn = self.outNCDir+"/"+str(var)+specific_date_string+".npy"
    #         np.save(fn, soil_water_state[var])
        
    #     crop_state = state['crop']
    #     for var in crop_state.iterItems():
    #         fn = self.outNCDir+"/"+str(var)+specific_date_string+".npy"
    #         np.save(fn, crop_state[var])
        
    # def resume(self):
    #     pass

    # def setState(self):
    #     logger.error("cannot set state")

    # def getState(self):
    #     result = {}
    #     result['groundwater'] = self.groundwater.getState()
    #     result['soilwater'] = self.soilwater.getState()
    #     result['landcover'] = self.landcover.getState()
    #     return result
    
    def getPseudoState(self):
        pass

    def getAllState(self):
        pass
