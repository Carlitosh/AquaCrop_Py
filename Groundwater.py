#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
from pcraster.framework import *
import pcraster as pcr

import logging
logger = logging.getLogger(__name__)

import VirtualOS as vos
# from ncConverter import *
# import ETPFunctions as refPotET

import logging
logger = logging.getLogger(__name__)

class Groundwater(object):

    def __init__(self, iniItems, landmask, spinUp):
        object.__init__(self)

        self.cloneMap = iniItems.cloneMap
        self.tmpDir = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']
        
        # landmask/area of interest
        self.landmask = landmask
        if iniItems.globalOptions['landmask'] != "None":
           self.landmask = vos.readPCRmapClone(iniItems.globalOptions['landmask'],\
                                               self.cloneMap,self.tmpDir,self.inputDir)

        # check if water table is to be modelled; if it is, check whether it is variable or constant
        # TODO: if a constant groundwater table is used then the groundwater 
        # NB added additional option to config: VariableWaterTable (= 1|0, i.e. True|False)

        self.WaterTable = bool(int(iniItems.groundwaterOptions['WaterTable']))
        self.VariableWaterTable = bool(int(iniItems.groundwaterOptions['VariableWaterTable']))

        if self.WaterTable:
            self.gwFileNC = iniItems.groundwaterOptions['groundwaterNC']
            self.gwVarName = iniItems.groundwaterOptions['groundwaterVariableName']

        # TODO: decide whether this is necessary
        # daily time step
        self.usingDailyTimeStepForcingData = False
        if iniItems.timeStep == 1.0 and iniItems.timeStepUnit == "day":
            self.usingDailyTimeStepForcingData = True

        # TODO: work out how to do this
        # option to use netcdf file that is defined per year (redundant?)
        self.groundwater_set_per_year  = False
            
        # make the iniItems available for the other modules:
        self.iniItems = iniItems
        
        # SM: ignore reporting for now (see meteo.py for example of how to do reporting)
        
    def update(self,currTimeStep):
        pass

    def read_forcings(self,currTimeStep):

        # TODO: decide whether this section is necessary
        # method for finding time indexes in the groundwater netdf file:
        # - the default one
        method_for_time_index = None
        # - based on the ini/configuration file (if given)

        if 'time_index_method_for_groundwater_netcdf' in self.iniItems.groundwaterOptions.keys() and\
                                                           self.iniItems.groundwaterOptions['time_index_method_for_groundwater_netcdf'] != "None":
            method_for_time_index = self.iniItems.groundwaterOptions['time_index_method_for_groundwater_netcdf']

        # reading groundwater:
        if self.WaterTable:
            
            if self.VariableWaterTable:
                if self.groundwater_set_per_year:
                    nc_file_per_year = self.gwFileNC %(float(currTimeStep.year), float(currTimeStep.year))
                    self.zGW = vos.netcdf2PCRobjClone(nc_file_per_year,\
                                                                     self.gwVarName,\
                                                                     str(currTimeStep.fulldate),\
                                                                     useDoy = method_for_time_index,\
                                                                     cloneMapFileName = self.cloneMap,\
                                                                     LatitudeLongitude = True)
                else:
                    self.zGW = vos.netcdf2PCRobjClone(self.gwFileNC,\
                                                                     self.gwVarName,\
                                                                     str(currTimeStep.fulldate),\
                                                                     useDoy = method_for_time_index,\
                                                                     cloneMapFileName = self.cloneMap,\
                                                                     LatitudeLongitude = True)
                    
            else:
                self.zGW = vos.netcdf2PCRobjCloneWithoutTime(self.gwFileNC,\
                                                                            self.gwVarName,\
                                                                            cloneMapFileName = self.cloneMap,\
                                                                            LatitudeLongitude = True)
                
        else:
            self.zGW = None  # TODO: check that this is OK in other parts of the program

        # SM: ignore reporting for now (see meteo.py for example of how to do reporting)
