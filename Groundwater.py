#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model
import os
import time
import numpy as np
import VirtualOS as vos

import logging
logger = logging.getLogger(__name__)

class Groundwater(object):

    def __init__(self, Groundwater_variable):
        self.var = Groundwater_variable

    def initial(self):
        self.var.WaterTable = bool(int(self.var._configuration.groundwaterOptions['WaterTable']))
        self.var.VariableWaterTable = bool(int(self.var._configuration.groundwaterOptions['VariableWaterTable']))
        self.var.DailyForcingData = bool(int(self.var._configuration.groundwaterOptions['DailyForcingData']))

        if self.var.WaterTable:
            self.var.gwFileNC = self.var._configuration.groundwaterOptions['groundwaterNC']
            self.var.gwVarName = self.var._configuration.groundwaterOptions['groundwaterVariableName']

        # # daily time step
        # self.var.usingDailyTimeStepForcingData = False
        # if self.var._configuration.timeStep == 1.0 and self.var._configuration.timeStepUnit == "day":
        #     self.var.usingDailyTimeStepForcingData = True

    def read(self):

        # method for finding time indexes in the groundwater netdf file:
        # - the default one
        method_for_time_index = None
        # - based on the ini/configuration file (if given)
        # if 'time_index_method_for_groundwater_netcdf' in self.var._configuration.groundwaterOptions.keys() and\
        #                                                    self.var._configuration.groundwaterOptions['time_index_method_for_groundwater_netcdf'] != "None":
        #     method_for_time_index = self.var._configuration.groundwaterOptions['time_index_method_for_groundwater_netcdf']
        
        # reading groundwater:
        if self.var.WaterTable:
            if self.var.VariableWaterTable:

                # DailyForcingData is a logical indicating whether a separate
                # netCDF is used for each time step - use this for coupling
                if self.var.DailyForcingData:

                    day, month, year = self.var._modelTime.day, self.var._modelTime.month, self.var._modelTime.year

                    # Fill named placeholders (NB we have already checked that
                    # the specified filename contains these placeholders)
                    gwFileNC = self.var.gwFileNC.format(day=day, month=month, year=year)

                    # Check whether the file is present in the filesystem; if
                    # it doesn't, enter a while loop which periodically checks
                    # whether the file exists. We specify a maximum wait time
                    # in order to prevent the model hanging if the file never
                    # materialises.
                    
                    exists = os.path.exists(gwFileNC)
                    max_wait_time = 60
                    wait_time = 0.1
                    total_wait_time = 0
                    while exists is False and total_wait_time <= max_wait_time:
                        time.sleep(wait_time)
                        exists = os.path.exists(gwFileNC)
                        total_wait_time += wait_time

                    if not exists:
                        print "groundwater file doesn't exist and maximum wait time exceeded"
                        raise   # TODO: make error class

                    self.var.zGW = vos.netcdf2PCRobjCloneWithoutTime(gwFileNC,
                                                                     self.var.gwVarName,
                                                                     cloneMapFileName = self.var.cloneMap,
                                                                     LatitudeLongitude = True)
                        
                else:
                    self.var.zGW = vos.netcdf2PCRobjClone(self.var.gwFileNC,
                                                          self.var.gwVarName,
                                                          str(currTimeStep.fulldate),
                                                          useDoy = method_for_time_index,
                                                          cloneMapFileName = self.var.cloneMap,
                                                          LatitudeLongitude = True)
                    
            else:
                self.var.zGW = vos.netcdf2PCRobjCloneWithoutTime(self.var.gwFileNC,
                                                                 self.var.gwVarName,
                                                                 cloneMapFileName = self.var.cloneMap,
                                                                 LatitudeLongitude = True)
                
        else:
            self.var.zGW = None
    
    def dynamic(self):
        self.read()
