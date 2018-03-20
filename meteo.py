#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
from pcraster.framework import *
import pcraster as pcr
import numpy as np

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos
# from ncConverter import *
# import ETPFunctions as refPotET

import logging
logger = logging.getLogger(__name__)

class Meteo(object):

    def __init__(self,iniItems,landmask,spinUp):
        object.__init__(self)

        self.cloneMap = iniItems.cloneMap
        self.tmpDir = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']
        
        # landmask/area of interest
        self.landmask = landmask
        if iniItems.globalOptions['landmask'] != "None":
           self.landmask = vos.readPCRmapClone(iniItems.globalOptions['landmask'],\
                                               self.cloneMap,self.tmpDir,self.inputDir)

        # # option to ignore snow (temperature will be set to 25 deg C if this option is activated)
        # self.ignore_snow = False
        # if 'ignoreSnow' in iniItems.meteoOptions.keys() and iniItems.meteoOptions['ignoreSnow'] == "True":
        #     self.ignore_snow = True

        self.preFileNC = iniItems.meteoOptions['precipitationNC']
        self.tmpFileNC = iniItems.meteoOptions['temperatureNC']
        self.etpFileNC = iniItems.meteoOptions['refETPotFileNC']

        # self.refETPotMethod = iniItems.meteoOptions['referenceETPotMethod']
        # if self.refETPotMethod == 'Hamon': self.latitudes = \
        #                                    pcr.ycoordinate(self.cloneMap) # needed to calculate 'referenceETPot'
        # if self.refETPotMethod == 'Input': self.etpFileNC = \
        #                      iniItems.meteoOptions['refETPotFileNC']              

        #-----------------------------------------------------------------------			
        # NOTE: RvB 13/07/2016 Added correction constant and factor and variable name
        # to allow for easier use of netCDF climate inpute files
        # EHS 20/08/2016 modified for more flexibilities.  
        # - meteo conversion factors
        self.preConst       = 0.0
        self.preFactor      = 1.0
        self.tmpConst       = 0.0
        self.tmpFactor      = 1.0
        self.refETPotConst  = 0.0
        self.refETPotFactor = 1.0
        # self.read_meteo_conversion_factors(iniItems.meteoOptions)

        # - variable names      
        self.preVarName      = 'precipitation' 
        self.tmnVarName      = 'Tmin'
        self.tmxVarName      = 'Tmax'
        self.refETPotVarName = 'evapotranspiration'

        # variable names
        self.read_meteo_variable_names(iniItems.meteoOptions)

        # daily time step
        self.usingDailyTimeStepForcingData = False
        if iniItems.timeStep == 1.0 and iniItems.timeStepUnit == "day":
            self.usingDailyTimeStepForcingData = True

        # SM: ignore downscaling for now
        # # forcing downscaling options:
        # self.forcingDownscalingOptions(iniItems)

        # option to use netcdf files that are defined per year (one file for each year)
        # self.precipitation_set_per_year  = iniItems.meteoOptions['precipitation_set_per_year'] == "True"
        # self.temperature_set_per_year    = iniItems.meteoOptions['temperature_set_per_year'] == "True"
        # self.refETPotFileNC_set_per_year = iniItems.meteoOptions['refETPotFileNC_set_per_year'] == "True" 
        self.precipitation_set_per_year  = False
        self.temperature_set_per_year    = False
        self.refETPotFileNC_set_per_year = False
        
        # make the iniItems available for the other modules:
        self.iniItems = iniItems

        # SM: ignore reporting for now
        
        # self.report = True
        # try:
        #     self.outDailyTotNC = iniItems.meteoOptions['outDailyTotNC'].split(",")
        #     self.outMonthTotNC = iniItems.meteoOptions['outMonthTotNC'].split(",")
        #     self.outMonthAvgNC = iniItems.meteoOptions['outMonthAvgNC'].split(",")
        #     self.outMonthEndNC = iniItems.meteoOptions['outMonthEndNC'].split(",")
        #     self.outAnnuaTotNC = iniItems.meteoOptions['outAnnuaTotNC'].split(",")
        #     self.outAnnuaAvgNC = iniItems.meteoOptions['outAnnuaAvgNC'].split(",")
        #     self.outAnnuaEndNC = iniItems.meteoOptions['outAnnuaEndNC'].split(",")
        # except:
        #     self.report = False
        # if self.report == True:
        #     # daily output in netCDF files:
        #     self.outNCDir  = iniItems.outNCDir
        #     self.netcdfObj = PCR2netCDF(iniItems)
        #     #
        #     if self.outDailyTotNC[0] != "None":
        #         for var in self.outDailyTotNC:
        #             # creating the netCDF files:
        #             self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
        #                                         str(var)+"_dailyTot.nc",\
        #                                             var,"undefined")
        #     # MONTHly output in netCDF files:
        #     # - cummulative
        #     if self.outMonthTotNC[0] != "None":
        #         for var in self.outMonthTotNC:
        #             # initiating monthlyVarTot (accumulator variable):
        #             vars(self)[var+'MonthTot'] = None
        #             # creating the netCDF files:
        #             self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
        #                                         str(var)+"_monthTot.nc",\
        #                                             var,"undefined")
        #     # - average
        #     if self.outMonthAvgNC[0] != "None":
        #         for var in self.outMonthAvgNC:
        #             # initiating monthlyTotAvg (accumulator variable)
        #             vars(self)[var+'MonthTot'] = None
        #             # initiating monthlyVarAvg:
        #             vars(self)[var+'MonthAvg'] = None
        #              # creating the netCDF files:
        #             self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
        #                                         str(var)+"_monthAvg.nc",\
        #                                             var,"undefined")
        #     # - last day of the month
        #     if self.outMonthEndNC[0] != "None":
        #         for var in self.outMonthEndNC:
        #              # creating the netCDF files:
        #             self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
        #                                         str(var)+"_monthEnd.nc",\
        #                                             var,"undefined")
        #     # YEARly output in netCDF files:
        #     # - cummulative
        #     if self.outAnnuaTotNC[0] != "None":
        #         for var in self.outAnnuaTotNC:
        #             # initiating yearly accumulator variable:
        #             vars(self)[var+'AnnuaTot'] = None
        #             # creating the netCDF files:
        #             self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
        #                                         str(var)+"_annuaTot.nc",\
        #                                             var,"undefined")
        #     # - average
        #     if self.outAnnuaAvgNC[0] != "None":
        #         for var in self.outAnnuaAvgNC:
        #             # initiating annualyVarAvg:
        #             vars(self)[var+'AnnuaAvg'] = None
        #             # initiating annualyTotAvg (accumulator variable)
        #             vars(self)[var+'AnnuaTot'] = None
        #              # creating the netCDF files:
        #             self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
        #                                         str(var)+"_annuaAvg.nc",\
        #                                             var,"undefined")
        #     # - last day of the year
        #     if self.outAnnuaEndNC[0] != "None":
        #         for var in self.outAnnuaEndNC:
        #              # creating the netCDF files:
        #             self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
        #                                         str(var)+"_annuaEnd.nc",\
        #                                             var,"undefined")


    # def read_meteo_conversion_factors(self, meteoOptions):

    #     # SM: PCRaster 'cover' command is used to fill in
    #     # missing values from one or more expressions (in
    #     # this case 0.0 or 1.0)
        
    #     if 'precipitationConstant' in meteoOptions: self.preConst       = pcr.cover(vos.readPCRmapClone(meteoOptions['precipitationConstant'], self.cloneMap, self.tmpDir, self.inputDir), 0.0)
    #     if 'precipitationFactor'   in meteoOptions: self.preFactor      = pcr.cover(vos.readPCRmapClone(meteoOptions['precipitationFactor'  ], self.cloneMap, self.tmpDir, self.inputDir), 1.0)
    #     if 'temperatureConstant'   in meteoOptions: self.tmpConst       = pcr.cover(vos.readPCRmapClone(meteoOptions['temperatureConstant'  ], self.cloneMap, self.tmpDir, self.inputDir), 0.0)
    #     if 'temperatureFactor'     in meteoOptions: self.tmpFactor      = pcr.cover(vos.readPCRmapClone(meteoOptions['temperatureFactor'    ], self.cloneMap, self.tmpDir, self.inputDir), 1.0)
    #     if 'referenceEPotConstant' in meteoOptions: self.refETPotConst  = pcr.cover(vos.readPCRmapClone(meteoOptions['referenceEPotConstant'], self.cloneMap, self.tmpDir, self.inputDir), 0.0)
    #     if 'referenceEPotFactor'   in meteoOptions: self.refETPotFactor = pcr.cover(vos.readPCRmapClone(meteoOptions['referenceEPotFactor'  ], self.cloneMap, self.tmpDir, self.inputDir), 1.0)


    def read_meteo_variable_names(self, meteoOptions):

        if 'precipitationVariableName' in meteoOptions: self.preVarName = meteoOptions['precipitationVariableName']
        # if 'temperatureVariableName'   in meteoOptions: self.tmpVarName      = meteoOptions['temperatureVariableName'  ]
        if 'tminVariableName'   in meteoOptions: self.tmnVarName        = meteoOptions['tminVariableName'  ]
        if 'tmaxVariableName'   in meteoOptions: self.tmxVarName        = meteoOptions['tmaxVariableName'  ]
        if 'refETPotVariableName' in meteoOptions: self.refETPotVarName = meteoOptions['refETPotVariableName']

    # SM: ignore downscaling for now - implement in future versions
    # def forcingDownscalingOptions(self, iniItems):

    #     self.downscalePrecipitationOption  = False
    #     self.downscaleTemperatureOption    = False
    #     self.downscaleReferenceETPotOption = False

    #     if 'meteoDownscalingOptions' in iniItems.allSections:

    #         # downscaling options
    #         if iniItems.meteoDownscalingOptions['downscalePrecipitation']  == "True":
    #             self.downscalePrecipitationOption  = True  
    #             logger.info("Precipitation forcing will be downscaled to the cloneMap resolution.")

    #         if iniItems.meteoDownscalingOptions['downscaleTemperature']    == "True":
    #             self.downscaleTemperatureOption    = True  
    #             logger.info("Temperature forcing will be downscaled to the cloneMap resolution.")

    #         if iniItems.meteoDownscalingOptions['downscaleReferenceETPot'] == "True" and self.refETPotMethod != 'Hamon':
    #             self.downscaleReferenceETPotOption = True 
    #             logger.info("Reference potential evaporation will be downscaled to the cloneMap resolution.")

    #         # Note that for the Hamon method: referencePotET will be calculated based on temperature,  
    #         # therefore, we do not have to downscale it (particularly if temperature is already provided at high resolution). 

    #     if self.downscalePrecipitationOption or\
    #        self.downscaleTemperatureOption   or\
    #        self.downscaleReferenceETPotOption:

    #         # creating anomaly DEM
    #         highResolutionDEM = vos.readPCRmapClone(\
    #            iniItems.meteoDownscalingOptions['highResolutionDEM'],
    #            self.cloneMap,self.tmpDir,self.inputDir)
    #         highResolutionDEM = pcr.cover(highResolutionDEM, 0.0)
    #         highResolutionDEM = pcr.max(highResolutionDEM, 0.0)
    #         self.meteoDownscaleIds = vos.readPCRmapClone(\
    #            iniItems.meteoDownscalingOptions['meteoDownscaleIds'],
    #            self.cloneMap,self.tmpDir,self.inputDir,isLddMap=False,cover=None,isNomMap=True)
    #         self.cellArea = vos.readPCRmapClone(\
    #            iniItems.routingOptions['cellAreaMap'],
    #            self.cloneMap,self.tmpDir,self.inputDir)
    #         loweResolutionDEM = pcr.areatotal(pcr.cover(highResolutionDEM*self.cellArea, 0.0),\
    #                                           self.meteoDownscaleIds)/\
    #                             pcr.areatotal(pcr.cover(self.cellArea, 0.0),\
    #                                           self.meteoDownscaleIds)                  
    #         self.anomalyDEM = highResolutionDEM - loweResolutionDEM    # unit: meter  

    #         # temperature lapse rate (netCDF) file 
    #         self.temperLapseRateNC = vos.getFullPath(iniItems.meteoDownscalingOptions[\
    #                                     'temperLapseRateNC'],self.inputDir)                         
    #         self.temperatCorrelNC  = vos.getFullPath(iniItems.meteoDownscalingOptions[\
    #                                     'temperatCorrelNC'],self.inputDir)                    # TODO: Remove this criteria.                         

    #         # precipitation lapse rate (netCDF) file 
    #         self.precipLapseRateNC = vos.getFullPath(iniItems.meteoDownscalingOptions[\
    #                                     'precipLapseRateNC'],self.inputDir)
    #         self.precipitCorrelNC  = vos.getFullPath(iniItems.meteoDownscalingOptions[\
    #                                     'precipitCorrelNC'],self.inputDir)                    # TODO: Remove this criteria.                           

    #     else:
    #         logger.info("No forcing downscaling is implemented.")

    #     # forcing smoothing options: - THIS is still experimental. PS: MUST BE TESTED.
    #     self.forcingSmoothing = False
    #     if 'meteoDownscalingOptions' in iniItems.allSections and \
    #        'smoothingWindowsLength' in iniItems.meteoDownscalingOptions.keys():

    #         if float(iniItems.meteoDownscalingOptions['smoothingWindowsLength']) > 0.0:
    #             self.forcingSmoothing = True
    #             self.smoothingWindowsLength = vos.readPCRmapClone(\
    #                iniItems.meteoDownscalingOptions['smoothingWindowsLength'],
    #                self.cloneMap,self.tmpDir,self.inputDir)
    #             msg = "Forcing data will be smoothed with 'windowaverage' using the window length:"+str(iniItems.meteoDownscalingOptions['smoothingWindowsLength'])
    #             logger.info(msg)   

    # SM: this function doesn't appear to be used anywhere
    # def perturb(self, name, **parameters):

    #     if name == "precipitation":

    #         # perturb the precipitation
    #         self.precipitation = self.precipitation * \
    #         pcr.min(pcr.max((1 + mapnormal() * parameters['standard_deviation']),0.01),2.0)
    #         #TODO: Please also make sure that precipitation >= 0
    #         #TODO: Add minimum and maximum 

    #     else:
    #         print("Error: only precipitation may be updated at this time")
    #         return -1


    def update(self,currTimeStep):
        #TODO: calculate  referencePotET
        pass

    # SM: as above, ignore downscaling options for now
    # def downscalePrecipitation(self, currTimeStep, useFactor = True, minCorrelationCriteria = 0.85):
        
    #     preSlope = 0.001 * vos.netcdf2PCRobjClone(\
    #                        self.precipLapseRateNC, 'precipitation',\
    #                        currTimeStep.month, useDoy = "Yes",\
    #                        cloneMapFileName=self.cloneMap,\
    #                        LatitudeLongitude = True)
    #     preSlope = pcr.cover(preSlope, 0.0)
    #     preSlope = pcr.max(0.,preSlope)
        
    #     preCriteria = vos.netcdf2PCRobjClone(\
    #                  self.precipitCorrelNC, 'precipitation',\
    #                  currTimeStep.month, useDoy = "Yes",\
    #                  cloneMapFileName=self.cloneMap,\
    #                  LatitudeLongitude = True)
    #     preSlope = pcr.ifthenelse(preCriteria > minCorrelationCriteria,\
    #                preSlope, 0.0)             
    #     preSlope = pcr.cover(preSlope, 0.0)
    
    #     if useFactor == True:
    #         factor = pcr.max(0.,self.precipitation + preSlope*self.anomalyDEM)
    #         factor = factor / \
    #                  pcr.areaaverage(factor, self.meteoDownscaleIds)
    #         factor = pcr.cover(factor, 1.0)
    #         self.precipitation = factor * self.precipitation
    #     else:
    #         self.precipitation = self.precipitation + preSlope*self.anomalyDEM

    #     self.precipitation = pcr.max(0.0, self.precipitation)

    # def downscaleTemperature(self, currTimeStep, useFactor = False, maxCorrelationCriteria = -0.75, zeroCelciusInKelvin = 273.15):
        
    #     tmpSlope = 1.000 * vos.netcdf2PCRobjClone(\
    #                        self.temperLapseRateNC, 'temperature',\
    #                        currTimeStep.month, useDoy = "Yes",\
    #                        cloneMapFileName=self.cloneMap,\
    #                        LatitudeLongitude = True)
    #     tmpSlope = pcr.min(0.,tmpSlope)  # must be negative
    #     tmpCriteria = vos.netcdf2PCRobjClone(\
    #                   self.temperatCorrelNC, 'temperature',\
    #                   currTimeStep.month, useDoy = "Yes",\
    #                   cloneMapFileName=self.cloneMap,\
    #                   LatitudeLongitude = True)
    #     tmpSlope = pcr.ifthenelse(tmpCriteria < maxCorrelationCriteria,\
    #                tmpSlope, 0.0)             
    #     tmpSlope = pcr.cover(tmpSlope, 0.0)
    
    #     if useFactor == True:
    #         temperatureInKelvin = self.temperature + zeroCelciusInKelvin
    #         factor = pcr.max(0.0, temperatureInKelvin + tmpSlope * self.anomalyDEM)
    #         factor = factor / \
    #                  pcr.areaaverage(factor, self.meteoDownscaleIds)
    #         factor = pcr.cover(factor, 1.0)
    #         self.temperature = factor * temperatureInKelvin - zeroCelciusInKelvin
    #     else:
    #         self.temperature = self.temperature + tmpSlope*self.anomalyDEM

    # def downscaleReferenceETPot(self, zeroCelciusInKelvin = 273.15):
        
    #     temperatureInKelvin = self.temperature + zeroCelciusInKelvin
    #     factor = pcr.max(0.0, temperatureInKelvin)
    #     factor = factor / \
    #              pcr.areaaverage(factor, self.meteoDownscaleIds)
    #     factor = pcr.cover(factor, 1.0)
    #     self.referencePotET = pcr.max(0.0, factor * self.referencePotET)

    def read_forcings(self,currTimeStep):

        #-----------------------------------------------------------------------
        # NOTE: RvB 13/07/2016 hard-coded reference to the variable names
        # preciptiation, temperature and evapotranspiration have been replaced
        # by the variable names used in the netCDF and passed from the ini file
        #-----------------------------------------------------------------------

        
        # method for finding time indexes in the precipitation netdf file:
        # - the default one
        method_for_time_index = None
        # - based on the ini/configuration file (if given)
        if 'time_index_method_for_precipitation_netcdf' in self.iniItems.meteoOptions.keys() and\
                                                           self.iniItems.meteoOptions['time_index_method_for_precipitation_netcdf'] != "None":
            method_for_time_index = self.iniItems.meteoOptions['time_index_method_for_precipitation_netcdf']
        
        # reading precipitation:
        if self.precipitation_set_per_year:
            #~ print currTimeStep.year
            nc_file_per_year = self.preFileNC %(float(currTimeStep.year), float(currTimeStep.year))
            self.precipitation = vos.netcdf2PCRobjClone(\
                                      nc_file_per_year, self.preVarName,\
                                      str(currTimeStep.fulldate), 
                                      useDoy = method_for_time_index,
                                      cloneMapFileName = self.cloneMap,\
                                      LatitudeLongitude = True)
        else:
            self.precipitation = vos.netcdf2PCRobjClone(self.preFileNC,\
                                                        self.preVarName,\
                                                        str(currTimeStep.fulldate),\
                                                        useDoy = method_for_time_index,\
                                                        cloneMapFileName = self.cloneMap,\
                                                        LatitudeLongitude = True)

        #-----------------------------------------------------------------------
        # NOTE: RvB 13/07/2016 added to automatically update precipitation

        # PCRaster.ifthen(condition, expression) - cell values of
        # condition (i.e. landmask, here) are interpreted as boolean
        # values where 1 = True and 0 = False...
        # (http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/manual/op_ifthen.html)

        # OLD:
        # self.precipitation  = self.preConst + self.preFactor * pcr.ifthen(self.landmask, self.precipitation)

        # NEW:
        # TODO: decided where np.nan is an appropriate missing value
        self.precipitation  = self.preConst + self.preFactor * np.where(self.landmask, self.precipitation, np.nan)

        #-----------------------------------------------------------------------

        # make sure that precipitation is always positive

        # OLD:
        # self.precipitation = pcr.max(0., self.precipitation)
        # self.precipitation = pcr.cover(  self.precipitation, 0.0)

        # NEW
        self.precipitation = np.maximum(0.0, self.precipitation)
        self.precipitation[np.isnan(self.precipitation)] = 0.0
        
        # ignore very small values of precipitation (less than 0.00001 m/day or less than 0.01 kg.m-2.day-1 )

        # OLD:
        # if self.usingDailyTimeStepForcingData:
        #     self.precipitation = pcr.rounddown(self.precipitation*100000.)/100000.

        # NEW:
        if self.usingDailyTimeStepForcingData:
            self.precipitation = np.floor(self.precipitation * 100000.)/100000.
        
        # method for finding time index in the temperature netdf file:
        # - the default one
        method_for_time_index = None
        # - based on the ini/configuration file (if given)
        if 'time_index_method_for_temperature_netcdf' in self.iniItems.meteoOptions.keys() and\
                                                         self.iniItems.meteoOptions['time_index_method_for_temperature_netcdf'] != "None":
            method_for_time_index = self.iniItems.meteoOptions['time_index_method_for_temperature_netcdf']

        # reading temperature
        # if self.temperature_set_per_year:
        #     nc_file_per_year = self.tmpFileNC %(int(currTimeStep.year), int(currTimeStep.year))
        #     self.temperature = vos.netcdf2PCRobjClone(\
        #                               nc_file_per_year, self.tmpVarName,\
        #                               str(currTimeStep.fulldate), 
        #                               useDoy = method_for_time_index,
        #                               cloneMapFileName = self.cloneMap,\
        #                               LatitudeLongitude = True)
        # else:
        #     self.temperature = vos.netcdf2PCRobjClone(\
        #                          self.tmpFileNC,self.tmpVarName,\
        #                          str(currTimeStep.fulldate), 
        #                          useDoy = method_for_time_index,
        #                          cloneMapFileName=self.cloneMap,\
        #                          LatitudeLongitude = True)
        if self.temperature_set_per_year:
            tmn_nc_file_per_year = self.tmpFileNC %(int(currTimeStep.year), int(currTimeStep.year))
            tmx_nc_file_per_year = self.tmpFileNC %(int(currTimeStep.year), int(currTimeStep.year))
            self.tmin = vos.netcdf2PCRobjClone(tmn_nc_file_per_year,\
                                               self.tmnVarName,\
                                               str(currTimeStep.fulldate),\
                                               useDoy = method_for_time_index,\
                                               cloneMapFileName = self.cloneMap,\
                                               LatitudeLongitude = True)

            self.tmax = vos.netcdf2PCRobjClone(tmx_nc_file_per_year,\
                                               self.tmxVarName,\
                                               str(currTimeStep.fulldate),\
                                               useDoy = method_for_time_index,\
                                               cloneMapFileName = self.cloneMap,\
                                               LatitudeLongitude = True)
        else:
            self.tmin = vos.netcdf2PCRobjClone(self.tmpFileNC,\
                                               self.tmnVarName,\
                                               str(currTimeStep.fulldate),\
                                               useDoy = method_for_time_index,\
                                               cloneMapFileName = self.cloneMap,\
                                               LatitudeLongitude = True)

            self.tmax = vos.netcdf2PCRobjClone(self.tmpFileNC,\
                                               self.tmxVarName,\
                                               str(currTimeStep.fulldate),\
                                               useDoy = method_for_time_index,\
                                               cloneMapFileName = self.cloneMap,\
                                               LatitudeLongitude = True)

        #-----------------------------------------------------------------------
        # NOTE: RvB 13/07/2016 added to automatically update temperature
        # self.temperature    = self.tmpConst + self.tmpFactor * pcr.ifthen(self.landmask, self.temperature)

        # OLD:
        # self.tmin    = self.tmpConst + self.tmpFactor * pcr.ifthen(self.landmask, self.tmin)
        # self.tmax    = self.tmpConst + self.tmpFactor * pcr.ifthen(self.landmask, self.tmax)

        # NEW:
        self.tmin    = self.tmpConst + self.tmpFactor * np.where(self.landmask, self.tmin, np.nan)
        self.tmax    = self.tmpConst + self.tmpFactor * np.where(self.landmask, self.tmax, np.nan)

        #-----------------------------------------------------------------------

        # SM: ignore downscaling for now
        # # Downscaling precipitation and temperature
        # if self.downscalePrecipitationOption: self.downscalePrecipitation(currTimeStep)
        # if self.downscaleTemperatureOption: self.downscaleTemperature(currTimeStep)

        # SM: only allow input
        # # calculate or obtain referencePotET
        # if self.refETPotMethod == 'Hamon': self.referencePotET = \
        #                           refPotET.HamonPotET(self.temperature,\
        #                                               currTimeStep.doy,\
        #                                               self.latitudes)
        # if self.refETPotMethod == 'Input': 

        # method for finding time indexes in the precipitation netdf file:
        # - the default one
        method_for_time_index = None
        # - based on the ini/configuration file (if given)
        if 'time_index_method_for_ref_pot_et_netcdf' in self.iniItems.meteoOptions.keys() and\
                                                        self.iniItems.meteoOptions['time_index_method_for_ref_pot_et_netcdf'] != "None":
            method_for_time_index = self.iniItems.meteoOptions['time_index_method_for_ref_pot_et_netcdf']

        if self.refETPotFileNC_set_per_year: 
            nc_file_per_year = self.etpFileNC %(int(currTimeStep.year), int(currTimeStep.year))
            self.referencePotET = vos.netcdf2PCRobjClone(\
                                  nc_file_per_year, self.refETPotVarName,\
                                  str(currTimeStep.fulldate), 
                                  useDoy = method_for_time_index,
                                  cloneMapFileName = self.cloneMap,\
                                  LatitudeLongitude = True)
        else:
            self.referencePotET = vos.netcdf2PCRobjClone(\
                                  self.etpFileNC,self.refETPotVarName,\
                                  str(currTimeStep.fulldate), 
                                  useDoy = method_for_time_index,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        #-----------------------------------------------------------------------
        # NOTE: RvB 13/07/2016 added to automatically update reference potential evapotranspiration

        # OLD:
        # self.referencePotET = self.refETPotConst + self.refETPotFactor * pcr.ifthen(self.landmask, self.referencePotET)

        # NEW:
        self.referencePotET = self.refETPotConst + self.refETPotFactor * np.where(self.landmask, self.referencePotET, np.nan)

        #-----------------------------------------------------------------------

        # SM: ignore downscaling for now
        # # Downscaling referenceETPot (based on temperature)
        # if self.downscaleReferenceETPotOption: self.downscaleReferenceETPot()

        # SM: ignore smoothing for now
        # # smoothing:
        # if self.forcingSmoothing == True:
        #     logger.debug("Forcing data are smoothed.")   
        #     self.precipitation  = pcr.windowaverage(self.precipitation , self.smoothingWindowsLength)
        #     self.temperature    = pcr.windowaverage(self.temperature   , self.smoothingWindowsLength)
        #     self.referencePotET = pcr.windowaverage(self.referencePotET, self.smoothingWindowsLength)
        
        # rounding temperature values to minimize numerical errors (note only to minimize, not remove)
        # self.temperature   = pcr.roundoff(self.temperature*1000.)/1000.

        # OLD:
        # self.tmin = pcr.roundoff(self.tmin * 1000.) / 1000.
        # self.tmax = pcr.roundoff(self.tmax * 1000.) / 1000.

        # NEW: (NB pcr and np behaviour is different - pcr: 0.5 -> 1.0, np: 0.5 -> 0)
        self.tmin = np.round(self.tmin * 1000.) / 1000.
        self.tmax = np.round(self.tmax * 1000.) / 1000.
        
        # SM: snow isn't relevant here (?)
        # # ignore snow by setting temperature to 25 deg C
        # if self.ignore_snow: self.temperature = pcr.spatial(pcr.scalar(25.))

        # define precipitation, temperature and referencePotET ONLY at landmask area (for reporting):

        # OLD:
        # self.precipitation  = pcr.ifthen(self.landmask, self.precipitation)
        # # self.temperature    = pcr.ifthen(self.landmask, self.temperature)
        # self.tmin = pcr.ifthen(self.landmask, self.tmin)
        # self.tmax = pcr.ifthen(self.landmask, self.tmax)
        # self.referencePotET = pcr.ifthen(self.landmask, self.referencePotET)

        # NEW:
        self.precipitation  = np.where(self.landmask, self.precipitation, np.nan)
        self.tmin           = np.where(self.landmask, self.tmin, np.nan)
        self.tmax           = np.where(self.landmask, self.tmax, np.nan)
        self.referencePotET = np.where(self.landmask, self.referencePotET, np.nan)
        
        # make sure precipitation and referencePotET are always positive:

        # OLD:
        # self.precipitation  = pcr.max(0.0, self.precipitation)
        # self.referencePotET = pcr.max(0.0, self.referencePotET)

        # NEW:
        self.precipitation  = np.maximum(0.0, self.precipitation)
        self.referencePotET = np.maximum(0.0, self.referencePotET)
        
        # SM: ignore reporting for now

        # if self.report == True:
        #     timeStamp = datetime.datetime(currTimeStep.year,\
        #                                   currTimeStep.month,\
        #                                   currTimeStep.day,\
        #                                   0)
        #     # writing daily output to netcdf files
        #     timestepPCR = currTimeStep.timeStepPCR
        #     if self.outDailyTotNC[0] != "None":
        #         for var in self.outDailyTotNC:
        #             self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
        #                                  str(var)+"_dailyTot.nc",\
        #                                  var,\
        #                   pcr2numpy(self.__getattribute__(var),vos.MV),\
        #                                  timeStamp,timestepPCR-1)

        #     # writing monthly output to netcdf files
        #     # -cummulative
        #     if self.outMonthTotNC[0] != "None":
        #         for var in self.outMonthTotNC:

        #             # introduce variables at the beginning of simulation or
        #             #     reset variables at the beginning of the month
        #             if currTimeStep.timeStepPCR == 1 or \
        #                currTimeStep.day == 1:\
        #                vars(self)[var+'MonthTot'] = pcr.scalar(0.0)

        #             # accumulating
        #             vars(self)[var+'MonthTot'] += vars(self)[var]

        #             # reporting at the end of the month:
        #             if currTimeStep.endMonth == True: 
        #                 self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
        #                                  str(var)+"_monthTot.nc",\
        #                                  var,\
        #                   pcr2numpy(self.__getattribute__(var+'MonthTot'),\
        #                    vos.MV),timeStamp,currTimeStep.monthIdx-1)
        #     # -average
        #     if self.outMonthAvgNC[0] != "None":
        #         for var in self.outMonthAvgNC:
        #             # only if a accumulator variable has not been defined: 
        #             if var not in self.outMonthTotNC: 

        #                 # introduce accumulator at the beginning of simulation or
        #                 #     reset accumulator at the beginning of the month
        #                 if currTimeStep.timeStepPCR == 1 or \
        #                    currTimeStep.day == 1:\
        #                    vars(self)[var+'MonthTot'] = pcr.scalar(0.0)
        #                 # accumulating
        #                 vars(self)[var+'MonthTot'] += vars(self)[var]

        #             # calculating average & reporting at the end of the month:
        #             if currTimeStep.endMonth == True:
        #                 vars(self)[var+'MonthAvg'] = vars(self)[var+'MonthTot']/\
        #                                              currTimeStep.day  
        #                 self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
        #                                  str(var)+"_monthAvg.nc",\
        #                                  var,\
        #                   pcr2numpy(self.__getattribute__(var+'MonthAvg'),\
        #                    vos.MV),timeStamp,currTimeStep.monthIdx-1)
        #     #
        #     # -last day of the month
        #     if self.outMonthEndNC[0] != "None":
        #         for var in self.outMonthEndNC:
        #             # reporting at the end of the month:
        #             if currTimeStep.endMonth == True: 
        #                 self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
        #                                  str(var)+"_monthEnd.nc",\
        #                                  var,\
        #                   pcr2numpy(self.__getattribute__(var),vos.MV),\
        #                                  timeStamp,currTimeStep.monthIdx-1)

        #     # writing yearly output to netcdf files
        #     # -cummulative
        #     if self.outAnnuaTotNC[0] != "None":
        #         for var in self.outAnnuaTotNC:

        #             # introduce variables at the beginning of simulation or
        #             #     reset variables at the beginning of the month
        #             if currTimeStep.timeStepPCR == 1 or \
        #                currTimeStep.doy == 1:\
        #                vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)

        #             # accumulating
        #             vars(self)[var+'AnnuaTot'] += vars(self)[var]

        #             # reporting at the end of the year:
        #             if currTimeStep.endYear == True: 
        #                 self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
        #                                  str(var)+"_annuaTot.nc",\
        #                                  var,\
        #                   pcr2numpy(self.__getattribute__(var+'AnnuaTot'),\
        #                    vos.MV),timeStamp,currTimeStep.annuaIdx-1)
        #     # -average
        #     if self.outAnnuaAvgNC[0] != "None":
        #         for var in self.outAnnuaAvgNC:
        #             # only if a accumulator variable has not been defined: 
        #             if var not in self.outAnnuaTotNC: 
        #                 # introduce accumulator at the beginning of simulation or
        #                 #     reset accumulator at the beginning of the year
        #                 if currTimeStep.timeStepPCR == 1 or \
        #                    currTimeStep.doy == 1:\
        #                    vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)
        #                 # accumulating
        #                 vars(self)[var+'AnnuaTot'] += vars(self)[var]
        #             #
        #             # calculating average & reporting at the end of the year:
        #             if currTimeStep.endYear == True:
        #                 vars(self)[var+'AnnuaAvg'] = vars(self)[var+'AnnuaTot']/\
        #                                              currTimeStep.doy  
        #                 self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
        #                                  str(var)+"_annuaAvg.nc",\
        #                                  var,\
        #                   pcr2numpy(self.__getattribute__(var+'AnnuaAvg'),\
        #                    vos.MV),timeStamp,currTimeStep.annuaIdx-1)
        #     #
        #     # -last day of the year
        #     if self.outAnnuaEndNC[0] != "None":
        #         for var in self.outAnnuaEndNC:
        #             # reporting at the end of the year:
        #             if currTimeStep.endYear == True: 
        #                 self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
        #                                  str(var)+"_annuaEnd.nc",\
        #                                  var,\
        #                   pcr2numpy(self.__getattribute__(var),vos.MV),\
        #                                  timeStamp,currTimeStep.annuaIdx-1)

        
