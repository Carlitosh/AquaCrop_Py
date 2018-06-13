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

    def initiate_reporting(self):
        """Function to create netCDF files for each output variable"""
        
        # output directory storing netcdf files:
        self.outNCDir  = self.configuration.outNCDir

        # object for reporting:
        self.netcdfObj = np2netCDF(self.configuration, self._model)

        # daily output in netCDF files:
        self.outDailyTotNC = ["None"]
        try:
            self.outDailyTotNC = list(set(self.configuration.reportingOptions['outDailyTotNC'].split(",")))
        except:
            pass
        
        if self.outDailyTotNC[0] != "None":
            for var in self.outDailyTotNC:
                
                logger.info("Creating the netcdf file for daily reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  

                dims       = varDicts.netcdf_dimensions[var]
                # includeDepthDimension = False
                # if 'depth' in dims:
                #     includeDepthDimension = True

                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+str(var)+"_dailyTot_output.nc",
                                            short_name,
                                            unit,
                                            dims,
                                            # includeDepthDimension,
                                            # includeRotationDimension,
                                            long_name)
        
        # list of variables that will be reported:
        self.variables_for_report = self.outDailyTotNC

    def post_processing(self):
        """Function to process model variables to output variables. In 
        most cases this simply involves copying model attributes
        """
        self.precipitation = self._model.meteo.precipitation
        self.th = self._model.soilwater.th
        self.Wr = self._model.soilwater.Wr
        self.zGW = self._model.groundwater.zGW
        self.SurfaceStorage = self._model.soilwater.SurfaceStorage
        self.Irr = self._model.soilwater.Irr
        self.Infl = self._model.soilwater.Infl
        self.Runoff = self._model.soilwater.Runoff
        self.DeepPerc = self._model.soilwater.DeepPerc
        self.CrTot = self._model.soilwater.CrTot
        self.GwIn = self._model.soilwater.GwIn
        self.EsAct = self._model.soilwater.EsAct
        self.Epot = self._model.soilwater.Epot
        self.TrAct = self._model.soilwater.TrAct
        self.Tpot = self._model.soilwater.Tpot
        self.GDD = self._model.landcover.GDD
        self.GDDcum = self._model.landcover.GDDcum
        self.Zroot = self._model.landcover.Zroot
        self.CC = self._model.landcover.CC
        self.CC_NS = self._model.landcover.CC_NS
        self.B = self._model.landcover.B
        self.B_NS = self._model.landcover.B_NS
        self.HI = self._model.landcover.HI
        self.HIadj = self._model.landcover.HIadj
        self.Y = self._model.landcover.Y

    def report(self):

        # recap all variables
        self.post_processing()

        # print np.max(self.th)
        
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
                dims       = varDicts.netcdf_dimensions[var]
                self.netcdfObj.data2NetCDF(self.outNCDir+"/"+str(var)+"_dailyTot_output.nc",
                                           short_name,
                                           dims,
                                           self.__getattribute__(var),
                                           # pcr.pcr2numpy(self.__getattribute__(var),vos.MV),
                                           timeStamp)
