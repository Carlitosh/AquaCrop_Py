#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
# from pcraster.framework import *
# import pcraster as pcr
import numpy as np

import VirtualOS as vos
# from ncConverter import *
# import ETPFunctions as refPotET

import logging
logger = logging.getLogger(__name__)

class Meteo(object):

    def __init__(self, Meteo_variable):
        self.var = Meteo_variable

    def initial(self):

        self.var.preFileNC = self.var._configuration.meteoOptions['precipitationNC']
        self.var.tmpFileNC = self.var._configuration.meteoOptions['temperatureNC']
        self.var.etpFileNC = self.var._configuration.meteoOptions['refETPotFileNC']

        # Meteo conversion factors
        self.var.preConst       = 0.0  # defaults
        self.var.preFactor      = 1.0
        self.var.tmpConst       = 0.0
        self.var.tmpFactor      = 1.0
        self.var.refETPotConst  = 0.0
        self.var.refETPotFactor = 1.0
        # TODO:
        # self.var.read_meteo_conversion_factors(_configuration.meteoOptions)

        # Variable names      
        self.var.preVarName      = 'precipitation'  # defaults
        self.var.tmnVarName      = 'Tmin'
        self.var.tmxVarName      = 'Tmax'
        self.var.refETPotVarName = 'evapotranspiration'
        self.read_meteo_variable_names(self.var._configuration.meteoOptions)

        # daily time step
        self.var.usingDailyTimeStepForcingData = False
        if self.var._configuration.timeStep == 1.0 and self.var._configuration.timeStepUnit == "day":
            self.var.usingDailyTimeStepForcingData = True

        # option to use netcdf files that are defined per year (one file for each year)
        # self.precipitation_set_per_year  = _configuration.meteoOptions['precipitation_set_per_year'] == "True"
        # self.temperature_set_per_year    = _configuration.meteoOptions['temperature_set_per_year'] == "True"
        # self.refETPotFileNC_set_per_year = _configuration.meteoOptions['refETPotFileNC_set_per_year'] == "True" 
        self.var.precipitation_set_per_year  = False
        self.var.temperature_set_per_year    = False
        self.var.refETPotFileNC_set_per_year = False

    def read_meteo_conversion_factors(self, meteoOptions):
        """Function to read conversion factors from configuration
        file
        """
        if 'precipitationConstant' in meteoOptions: self.var.preConst = meteoOptions['precipitationConstant']
        if 'precipitationFactor' in meteoOptions: self.var.preFactor = meteoOptions['precipitationFactor']
        if 'temperatureConstant' in meteoOptions: self.var.tmpConst = meteoOptions['temperatureConstant']
        if 'temperatureFactor' in meteoOptions: self.var.tmpFactor = meteoOptions['temperatureFactor']
        if 'ETpotConstant' in meteoOptions: self.var.refETPotConst = meteoOptions['ETpotConstant']
        if 'ETpotFactor' in meteoOptions: self.var.refETPotFactor = meteoOptions['ETpotFactor']
        
    def read_meteo_variable_names(self, meteoOptions):
        """Function to read netCDF variable names from the 
        configuration file
        """
        if 'precipitationVariableName' in meteoOptions: self.preVarName = meteoOptions['precipitationVariableName']
        if 'tminVariableName' in meteoOptions: self.tmnVarName = meteoOptions['tminVariableName'  ]
        if 'tmaxVariableName' in meteoOptions: self.tmxVarName = meteoOptions['tmaxVariableName'  ]
        if 'refETPotVariableName' in meteoOptions: self.refETPotVarName = meteoOptions['refETPotVariableName']

    def dynamic(self):
        
        # method for finding time indexes in the precipitation netdf file:
        # - the default one
        method_for_time_index = None
        # - based on the ini/configuration file (if given)
        if 'time_index_method_for_precipitation_netcdf' in self.var._configuration.meteoOptions.keys() and self.var._configuration.meteoOptions['time_index_method_for_precipitation_netcdf'] != "None":
            method_for_time_index = self.var._configuration.meteoOptions['time_index_method_for_precipitation_netcdf']
        
        # reading precipitation:
        if self.var.precipitation_set_per_year:
            nc_file_per_year = self.var.preFileNC %(float(self.var._modelTime.year), float(self.var._modelTime.year))
            self.var.precipitation = vos.netcdf2PCRobjClone(\
                                      nc_file_per_year, self.var.preVarName,\
                                      str(self.var._modelTime.fulldate), 
                                      useDoy = method_for_time_index,
                                      cloneMapFileName = self.var.cloneMap,\
                                      LatitudeLongitude = True)
        else:
            self.var.precipitation = vos.netcdf2PCRobjClone(self.var.preFileNC,\
                                                            self.var.preVarName,\
                                                            str(self.var._modelTime.fulldate),\
                                                            useDoy = method_for_time_index,\
                                                            cloneMapFileName = self.var.cloneMap,\
                                                            LatitudeLongitude = True)

        # TODO: decided where np.nan is an appropriate missing value
        self.var.precipitation  = self.var.preConst + self.var.preFactor * np.where(self.var.landmask, self.var.precipitation, np.nan)

        # make sure that precipitation is always positive
        self.var.precipitation = np.maximum(0.0, self.var.precipitation)
        self.var.precipitation[np.isnan(self.var.precipitation)] = 0.0
        
        # ignore very small values of precipitation (less than 0.00001 m/day or less than 0.01 kg.m-2.day-1 )
        if self.var.usingDailyTimeStepForcingData:
            self.var.precipitation = np.floor(self.var.precipitation * 100000.)/100000.
        
        # method for finding time index in the temperature netdf file:
        # - the default one
        method_for_time_index = None
        # - based on the ini/configuration file (if given)
        if 'time_index_method_for_temperature_netcdf' in self.var._configuration.meteoOptions.keys() and\
                                                         self.var._configuration.meteoOptions['time_index_method_for_temperature_netcdf'] != "None":
            method_for_time_index = self.var._configuration.meteoOptions['time_index_method_for_temperature_netcdf']

        # reading temperature
        if self.var.temperature_set_per_year:
            tmn_nc_file_per_year = self.var.tmpFileNC %(int(self.var._modelTime.year), int(self.var._modelTime.year))
            tmx_nc_file_per_year = self.var.tmpFileNC %(int(self.var._modelTime.year), int(self.var._modelTime.year))
            self.var.tmin = vos.netcdf2PCRobjClone(tmn_nc_file_per_year,\
                                               self.var.tmnVarName,\
                                               str(self.var._modelTime.fulldate),\
                                               useDoy = method_for_time_index,\
                                               cloneMapFileName = self.var.cloneMap,\
                                               LatitudeLongitude = True)

            self.var.tmax = vos.netcdf2PCRobjClone(tmx_nc_file_per_year,\
                                               self.var.tmxVarName,\
                                               str(self.var._modelTime.fulldate),\
                                               useDoy = method_for_time_index,\
                                               cloneMapFileName = self.var.cloneMap,\
                                               LatitudeLongitude = True)
        else:
            self.var.tmin = vos.netcdf2PCRobjClone(self.var.tmpFileNC,\
                                               self.var.tmnVarName,\
                                               str(self.var._modelTime.fulldate),\
                                               useDoy = method_for_time_index,\
                                               cloneMapFileName = self.var.cloneMap,\
                                               LatitudeLongitude = True)

            self.var.tmax = vos.netcdf2PCRobjClone(self.var.tmpFileNC,\
                                               self.var.tmxVarName,\
                                               str(self.var._modelTime.fulldate),\
                                               useDoy = method_for_time_index,\
                                               cloneMapFileName = self.var.cloneMap,\
                                               LatitudeLongitude = True)

        self.var.tmin = self.var.tmpConst + self.var.tmpFactor * np.where(self.var.landmask, self.var.tmin, np.nan)
        self.var.tmax = self.var.tmpConst + self.var.tmpFactor * np.where(self.var.landmask, self.var.tmax, np.nan)

        # round to nearest mm
        self.var.tmin = np.round(self.var.tmin * 1000.) / 1000.
        self.var.tmax = np.round(self.var.tmax * 1000.) / 1000.

        if 'time_index_method_for_ref_pot_et_netcdf' in self.var._configuration.meteoOptions.keys() and self.var._configuration.meteoOptions['time_index_method_for_ref_pot_et_netcdf'] != "None":
            method_for_time_index = self.var._configuration.meteoOptions['time_index_method_for_ref_pot_et_netcdf']

        if self.var.refETPotFileNC_set_per_year: 
            nc_file_per_year = self.var.etpFileNC %(int(self.var._modelTime.year), int(self.var._modelTime.year))
            self.var.referencePotET = vos.netcdf2PCRobjClone(
                nc_file_per_year,
                self.var.refETPotVarName,
                str(self.var._modelTime.fulldate), 
                useDoy = method_for_time_index,
                cloneMapFileName = self.var.cloneMap,
                LatitudeLongitude = True)
        else:
            self.var.referencePotET = vos.netcdf2PCRobjClone(
                self.var.etpFileNC,
                self.var.refETPotVarName,
                str(self.var._modelTime.fulldate), 
                useDoy = method_for_time_index,
                cloneMapFileName=self.var.cloneMap,
                LatitudeLongitude = True)

        self.var.referencePotET = self.var.refETPotConst + self.var.refETPotFactor * np.where(self.var.landmask, self.var.referencePotET, np.nan)
