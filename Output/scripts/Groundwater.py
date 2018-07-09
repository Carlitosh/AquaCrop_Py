#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np
import VirtualOS as vos

from crop_growth_funs import *

import logging
logger = logging.getLogger(__name__)

class Groundwater(object):

    def __init__(self, Groundwater_variable):
        self.var = Groundwater_variable

    def initial(self):

        # NB added additional option to config: VariableWaterTable (= 1|0, i.e. True|False)

        self.var.WaterTable = bool(int(self.var._configuration.groundwaterOptions['WaterTable']))
        self.var.VariableWaterTable = bool(int(self.var._configuration.groundwaterOptions['VariableWaterTable']))

        if self.var.WaterTable:
            self.var.gwFileNC = self.var._configuration.groundwaterOptions['groundwaterNC']
            self.var.gwVarName = self.var._configuration.groundwaterOptions['groundwaterVariableName']

        # daily time step
        self.var.usingDailyTimeStepForcingData = False
        if self.var._configuration.timeStep == 1.0 and self.var._configuration.timeStepUnit == "day":
            self.var.usingDailyTimeStepForcingData = True

        # option to use netcdf file that is defined per year (redundant?)
        self.var.groundwater_set_per_year  = False
            
    # def update(self,currTimeStep):
    #     pass

    def dynamic(self):

        # TODO: write coupling method for AMBHAS
        
        # method for finding time indexes in the groundwater netdf file:
        # - the default one
        method_for_time_index = None
        # - based on the ini/configuration file (if given)

        if 'time_index_method_for_groundwater_netcdf' in self.var._configuration.groundwaterOptions.keys() and\
                                                           self.var._configuration.groundwaterOptions['time_index_method_for_groundwater_netcdf'] != "None":
            method_for_time_index = self.var._configuration.groundwaterOptions['time_index_method_for_groundwater_netcdf']

        # reading groundwater:
        if self.var.WaterTable:
            
            if self.var.VariableWaterTable:
                if self.var.groundwater_set_per_year:
                    nc_file_per_year = self.var.gwFileNC %(float(currTimeStep.year), float(currTimeStep.year))
                    self.var.zGW = vos.netcdf2PCRobjClone(nc_file_per_year,\
                                                      self.var.gwVarName,\
                                                      str(currTimeStep.fulldate),\
                                                      useDoy = method_for_time_index,\
                                                      cloneMapFileName = self.var.cloneMap,\
                                                      LatitudeLongitude = True)
                else:
                    self.var.zGW = vos.netcdf2PCRobjClone(self.var.gwFileNC,\
                                                      self.var.gwVarName,\
                                                      str(currTimeStep.fulldate),\
                                                      useDoy = method_for_time_index,\
                                                      cloneMapFileName = self.var.cloneMap,\
                                                      LatitudeLongitude = True)
                    
            else:
                self.var.zGW = vos.netcdf2PCRobjCloneWithoutTime(self.var.gwFileNC,\
                                                             self.var.gwVarName,\
                                                             cloneMapFileName = self.var.cloneMap,\
                                                             LatitudeLongitude = True)
                
        else:
            self.var.zGW = None # TODO: check that this is OK in other parts of the program
