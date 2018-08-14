#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil

from ncConverter import *
from types import NoneType
import variable_list as varDicts

import logging
logger = logging.getLogger(__name__)

class Reporting(object):

    def __init__(self, configuration, model, modelTime):

        # model (e.g. PCR-GLOBWB) and modelTime object
        self._model = model
        self._modelTime = modelTime

        # configuration/setting from the ini file
        self.configuration = configuration
        
        # initiate reporting tool/object and its configuration
        self.initiate_reporting()

    def create_netcdf_file(self, var, suffix):
        
        short_name = varDicts.netcdf_short_name[var]
        unit       = varDicts.netcdf_unit[var]      
        long_name  = varDicts.netcdf_long_name[var]
        if long_name == None: long_name = short_name  
        dims       = varDicts.netcdf_dimensions[var]

        # creating netCDF files:
        self.netcdfObj.createNetCDF(self.outNCDir+"/"+str(var)+str(suffix)+".nc",
                                    short_name,
                                    unit,
                                    dims,
                                    long_name)

    def initiate_reporting(self):
        """Function to create netCDF files for each output variable"""
        
        # output directory storing netcdf files:
        self.outNCDir  = self.configuration.outNCDir

        # object for reporting:
        self.netcdfObj = np2netCDF(self.configuration, self._model)

        # daily output in netCDF files:
        # #############################        
        self.outDailyTotNC = ["None"]
        try:
            self.outDailyTotNC = list(set(self.configuration.reportingOptions['outDailyTotNC'].split(",")))
        except:
            pass
        
        if self.outDailyTotNC[0] != "None":
            for var in self.outDailyTotNC:                
                logger.info("Creating the netcdf file for reporting the daily value of variable %s.", str(var))
                self.create_netcdf_file(var, "_dailyTot_output")

        # month average in netCDF files:
        # ##############################
        self.outMonthAvgNC = ["None"]
        try:
            self.outMonthAvgNC = list(set(self.configuration.reportingOptions['outMonthAvgNC'].split(",")))
        except:
            pass
        if self.outMonthAvgNC[0] != "None":
            for var in self.outMonthAvgNC:
                logger.info("Creating the netcdf file for reporting the monthly average of variable %s.", str(var))
                self.create_netcdf_file(var, "_monthAvg_output")
                vars(self)[var+'_monthAvg'] = np.zeros((vars(self._model)[var].shape))
                
        # month end in netCDF files:
        # ##########################
        self.outMonthEndNC = ["None"]
        try:
            self.outMonthEndNC = list(set(self.configuration.reportingOptions['outMonthEndNC'].split(",")))
        except:
            pass
        if self.outMonthEndNC[0] != "None":
            for var in self.outMonthEndNC:
                logger.info("Creating the netcdf file for reporting the month end value of variable %s.", str(var))
                self.create_netcdf_file(var, "_monthEnd_output")
                vars(self)[var+'_monthEnd'] = np.zeros((vars(self._model)[var].shape))

        # month total in netCDF files:
        # ############################
        self.outMonthTotNC = ["None"]
        try:
            self.outMonthTotNC = list(set(self.configuration.reportingOptions['outMonthTotNC'].split(",")))
        except:
            pass
        if self.outMonthTotNC[0] != "None":
            for var in self.outMonthTotNC:
                logger.info("Creating the netcdf file for reporting the monthly average of variable %s.", str(var))
                self.create_netcdf_file(var, "_monthTot_output")
                vars(self)[var+'_monthTot'] = np.zeros((vars(self._model)[var].shape))

        # month maximum in netCDF files:
        # ##############################
        self.outMonthMaxNC = ["None"]
        try:
            self.outMonthMaxNC = list(set(self.configuration.reportingOptions['outMonthMaxNC'].split(",")))
        except:
            pass
        if self.outMonthMaxNC[0] != "None":
            for var in self.outMonthMaxNC:
                logger.info("Creating the netcdf file for reporting the monthly maximum of variable %s.", str(var))
                self.create_netcdf_file(var, "_monthMax_output")
                vars(self)[var+'_monthMax'] = np.zeros((vars(self._model)[var].shape))

        # year average in netCDF files:
        # ##############################
        self.outYearAvgNC = ["None"]
        try:
            self.outYearAvgNC = list(set(self.configuration.reportingOptions['outYearAvgNC'].split(",")))
        except:
            pass
        if self.outYearAvgNC[0] != "None":
            for var in self.outYearAvgNC:
                logger.info("Creating the netcdf file for reporting the yearly average of variable %s.", str(var))
                self.create_netcdf_file(var, "_yearAvg_output")
                vars(self)[var+'_yearAvg'] = np.zeros((vars(self._model)[var].shape))
                
        # year end in netCDF files:
        # ##########################
        self.outYearEndNC = ["None"]
        try:
            self.outYearEndNC = list(set(self.configuration.reportingOptions['outYearEndNC'].split(",")))
        except:
            pass
        if self.outYearEndNC[0] != "None":
            for var in self.outYearEndNC:
                logger.info("Creating the netcdf file for reporting the year end value of variable %s.", str(var))
                self.create_netcdf_file(var, "_yearEnd_output")
                vars(self)[var+'_yearEnd'] = np.zeros((vars(self._model)[var].shape))

        # year total in netCDF files:
        # ############################
        self.outYearTotNC = ["None"]
        try:
            self.outYearTotNC = list(set(self.configuration.reportingOptions['outYearTotNC'].split(",")))
        except:
            pass
        if self.outYearTotNC[0] != "None":
            for var in self.outYearTotNC:
                logger.info("Creating the netcdf file for reporting the yearly total of variable %s.", str(var))
                self.create_netcdf_file(var, "_yearTot_output")
                vars(self)[var+'_yearTot'] = np.zeros((vars(self._model)[var].shape))

        # year maximum in netCDF files:
        # ##############################
        self.outYearMaxNC = ["None"]
        try:
            self.outYearMaxNC = list(set(self.configuration.reportingOptions['outYearMaxNC'].split(",")))
        except:
            pass
        if self.outYearMaxNC[0] != "None":
            for var in self.outYearMaxNC:
                logger.info("Creating the netcdf file for reporting the yearly maximum of variable %s.", str(var))
                self.create_netcdf_file(var, "_yearMax_output")
                vars(self)[var+'_yearMax'] = np.zeros((vars(self._model)[var].shape))
                
        # list of variables that will be reported:
        self.variables_for_report = (
            self.outDailyTotNC +
            self.outMonthAvgNC +
            self.outMonthEndNC +
            self.outMonthTotNC +
            self.outMonthMaxNC +
            self.outYearAvgNC +
            self.outYearEndNC +
            self.outYearTotNC +
            self.outYearMaxNC
            )

        # reduce above list to unique values, and remove None
        self.variables_for_report = list(set(self.variables_for_report))                                             
        if "None" in self.variables_for_report: self.variables_for_report.remove("None")

    def post_processing(self):
        """Function to process model variables to output variables. In 
        most cases this simply involves copying model attributes
        """
        # if self.outDailyTotNC[0] != "None":
        #     for var in self.outDailyTotNC:
        #         vars(self)[var] = vars(self._model)[var]
        if len(self.variables_for_report) > 0:
            for var in self.variables_for_report:
                vars(self)[var] = vars(self._model)[var]
               
    def report(self):

        # recap all variables
        self.post_processing()

        # time stamp for reporting
        timeStamp = datetime.datetime(self._modelTime.year,\
                                      self._modelTime.month,\
                                      self._modelTime.day,\
                                      0)

        logger.info("reporting for time %s", self._modelTime.currTime)

        # writing daily output to netcdf files
        if self.outDailyTotNC[0] != "None":
            for var in self.outDailyTotNC:
                short_name = varDicts.netcdf_short_name[var]
                dims = varDicts.netcdf_dimensions[var]
                self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_dailyTot_output.nc",
                                           short_name,
                                           dims,
                                           self.__getattribute__(var),
                                           timeStamp)

        if self.outMonthAvgNC[0] != "None":
            for var in self.outMonthAvgNC:
                if self._modelTime.day == 1: vars(self)[var+'_monthAvg'].fill(0)
                vars(self)[var+'_monthAvg'] += vars(self)[var]
                if self._modelTime.endMonth:
                    vars(self)[var+'_monthAvg'] = vars(self)[var+'_monthAvg'] / self._modelTime.day
                    short_name = varDicts.netcdf_short_name[var]
                    dims = varDicts.netcdf_dimensions[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_monthAvg_output.nc",
                                               short_name,
                                               dims,
                                               self.__getattribute__(var+'_monthAvg'),
                                               timeStamp)
                    

        if self.outMonthEndNC[0] != "None":
            for var in self.outMonthEndNC:
                if self._modelTime.endMonth:
                    vars(self)[var+'_monthEnd'] = vars(self)[var]
                    short_name = varDicts.netcdf_short_name[var]
                    dims = varDicts.netcdf_dimensions[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_monthEnd_output.nc",
                                               short_name,
                                               dims,
                                               self.__getattribute__(var+'_monthEnd'),
                                               timeStamp)
                    
        if self.outMonthTotNC[0] != "None":
            for var in self.outMonthTotNC:
                if self._modelTime.day == 1: vars(self)[var+'_monthAvg'].fill(0)
                vars(self)[var+'_monthAvg'] += vars(self)[var]
                if self._modelTime.endMonth:
                    short_name = varDicts.netcdf_short_name[var]
                    dims = varDicts.netcdf_dimensions[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_monthTot_output.nc",
                                               short_name,
                                               dims,
                                               self.__getattribute__(var+'_monthAvg'),
                                               timeStamp)                    

        if self.outMonthMaxNC[0] != "None":
            for var in self.outMonthMaxNC:
                if self._modelTime.day == 1: vars(self)[var+'_monthMax'].fill(0)
                vars(self)[var+'_monthMax'].clip(vars(self)[var], None)
                if self._modelTime.endMonth:
                    short_name = varDicts.netcdf_short_name[var]
                    dims = varDicts.netcdf_dimensions[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_monthMax_output.nc",
                                               short_name,
                                               dims,
                                               self.__getattribute__(var+'_monthMax'),
                                               timeStamp)                    
                
        if self.outYearAvgNC[0] != "None":
            for var in self.outYearAvgNC:
                if self._modelTime.doy == 1: vars(self)[var+'_yearAvg'].fill(0)
                vars(self)[var+'_yearAvg'] += vars(self)[var]
                if self._modelTime.endYear:
                    vars(self)[var+'_yearAvg'] = vars(self)[var+'_yearAvg'] / self._modelTime.day
                    short_name = varDicts.netcdf_short_name[var]
                    dims = varDicts.netcdf_dimensions[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_yearAvg_output.nc",
                                               short_name,
                                               dims,
                                               self.__getattribute__(var+'_yearAvg'),
                                               timeStamp)
                    

        if self.outYearEndNC[0] != "None":
            for var in self.outYearEndNC:
                if self._modelTime.endYear:
                    vars(self)[var+'_yearEnd'] = vars(self)[var]
                    short_name = varDicts.netcdf_short_name[var]
                    dims = varDicts.netcdf_dimensions[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_yearEnd_output.nc",
                                               short_name,
                                               dims,
                                               self.__getattribute__(var+'_yearEnd'),
                                               timeStamp)
                    
        if self.outYearTotNC[0] != "None":
            for var in self.outYearTotNC:
                if self._modelTime.doy == 1: vars(self)[var+'_yearAvg'].fill(0)
                vars(self)[var+'_yearAvg'] += vars(self)[var]
                if self._modelTime.endYear:
                    short_name = varDicts.netcdf_short_name[var]
                    dims = varDicts.netcdf_dimensions[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_yearTot_output.nc",
                                               short_name,
                                               dims,
                                               self.__getattribute__(var+'_yearAvg'),
                                               timeStamp)                    

        if self.outYearMaxNC[0] != "None":
            for var in self.outYearMaxNC:
                if self._modelTime.doy == 1: vars(self)[var+'_yearMax'].fill(0)
                vars(self)[var+'_yearMax'] = np.clip(vars(self)[var+'_yearMax'], vars(self)[var], None)
                if self._modelTime.endYear:
                    short_name = varDicts.netcdf_short_name[var]
                    dims = varDicts.netcdf_dimensions[var]
                    # print self.__getattribute__(var).shape
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_yearMax_output.nc",
                                               short_name,
                                               dims,
                                               # self.__getattribute__(var),
                                               self.__getattribute__(var+'_yearMax'),
                                               timeStamp)                    
            
