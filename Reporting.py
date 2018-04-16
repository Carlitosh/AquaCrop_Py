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

        # # option for debugging to PCR-GLOBWB version 1.0
        # self.debug_to_version_one = False
        # if self.configuration.debug_to_version_one: self.debug_to_version_one = True

    def initiate_reporting(self):
        
        # output directory storing netcdf files:
        self.outNCDir  = self.configuration.outNCDir

        # object for reporting:
        self.netcdfObj = PCR2netCDF(self.configuration)

        # initiating netcdf files for reporting

        # - daily output in netCDF files:
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
                
                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_dailyTot_output.nc",\
                                            short_name,unit,long_name)

        # # - Monthly output in netCDF files (cummulative)
        # self.outMonthTotNC = ["None"]
        # try:
        #     self.outMonthTotNC = list(set(self.configuration.reportingOptions['outMonthTotNC'].split(",")))
        # except:
        #     pass
        # if self.outMonthTotNC[0] != "None":
        #     for var in self.outMonthTotNC:

        #         # initiating monthlyVarTot (accumulator variable):
        #         vars(self)[var+'MonthTot'] = None

        #         logger.info("Creating the netcdf file for monthly accumulation reporting for variable %s.", str(var))

        #         short_name = varDicts.netcdf_short_name[var]
        #         unit       = varDicts.netcdf_monthly_total_unit[var]      
        #         long_name  = varDicts.netcdf_long_name[var]
        #         if long_name == None: long_name = short_name  

        #         # creating netCDF files:
        #         self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                     "_monthTot_output.nc",\
        #                                     short_name,unit,long_name)

        # # -- average
        # self.outMonthAvgNC = ["None"]
        # try:
        #     self.outMonthAvgNC = list(set(self.configuration.reportingOptions['outMonthAvgNC'].split(",")))
        # except:
        #     pass
        # if self.outMonthAvgNC[0] != "None":

        #     for var in self.outMonthAvgNC:

        #         # initiating monthlyTotAvg (accumulator variable)
        #         vars(self)[var+'MonthTot'] = None

        #         # initiating monthlyVarAvg:
        #         vars(self)[var+'MonthAvg'] = None

        #         logger.info("Creating the netcdf file for monthly average reporting for variable %s.", str(var))

        #         short_name = varDicts.netcdf_short_name[var]
        #         unit       = varDicts.netcdf_unit[var]      
        #         long_name  = varDicts.netcdf_long_name[var]
        #         if long_name == None: long_name = short_name  

        #         # creating netCDF files:
        #         self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                     "_monthAvg_output.nc",\
        #                                     short_name,unit,long_name)
        # #
        # # -- last day of the month
        # self.outMonthEndNC = ["None"]
        # try:
        #     self.outMonthEndNC = list(set(self.configuration.reportingOptions['outMonthEndNC'].split(",")))
        # except:
        #     pass
        # if self.outMonthEndNC[0] != "None":

        #     for var in self.outMonthEndNC:

        #         logger.info("Creating the netcdf file for monthly end reporting for variable %s.", str(var))

        #         short_name = varDicts.netcdf_short_name[var]
        #         unit       = varDicts.netcdf_unit[var]      
        #         long_name  = varDicts.netcdf_long_name[var]
        #         if long_name == None: long_name = short_name  

        #         # creating netCDF files:
        #         self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                     "_monthEnd_output.nc",\
        #                                     short_name,unit,long_name)
        # #
        # # -- maximum of the month
        # self.outMonthMaxNC = ["None"]
        # try:
        #     self.outMonthMaxNC = list(set(self.configuration.reportingOptions['outMonthMaxNC'].split(",")))
        # except:
        #     pass
        # if self.outMonthMaxNC[0] != "None":

        #     for var in self.outMonthMaxNC:

        #         logger.info("Creating the netcdf file for monthly maximum reporting for variable %s.", str(var))

        #         short_name = varDicts.netcdf_short_name[var]
        #         unit       = varDicts.netcdf_unit[var]      
        #         long_name  = varDicts.netcdf_long_name[var]
        #         if long_name == None: long_name = short_name  

        #         # creating netCDF files:
        #         self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                     "_monthMax_output.nc",\
        #                                     short_name,unit,long_name)

        # # - YEARly output in netCDF files:
        # # -- cummulative
        # self.outAnnuaTotNC = ["None"]
        # try:
        #     self.outAnnuaTotNC = list(set(self.configuration.reportingOptions['outAnnuaTotNC'].split(",")))
        # except:
        #     pass
        # if self.outAnnuaTotNC[0] != "None":

        #     for var in self.outAnnuaTotNC:

        #         # initiating yearly accumulator variable:
        #         vars(self)[var+'AnnuaTot'] = None

        #         logger.info("Creating the netcdf file for annual accumulation reporting for variable %s.", str(var))

        #         short_name = varDicts.netcdf_short_name[var]
        #         unit       = varDicts.netcdf_yearly_total_unit[var]      
        #         long_name  = varDicts.netcdf_long_name[var]
        #         if long_name == None: long_name = short_name  

        #         # creating netCDF files:
        #         self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                     "_annuaTot_output.nc",\
        #                                     short_name,unit,long_name)
        # #
        # # -- average
        # self.outAnnuaAvgNC = ["None"]
        # try:
        #     self.outAnnuaAvgNC = list(set(self.configuration.reportingOptions['outAnnuaAvgNC'].split(",")))
        # except:
        #     pass
        # if self.outAnnuaAvgNC[0] != "None":

        #     for var in self.outAnnuaAvgNC:

        #         # initiating annualyVarAvg:
        #         vars(self)[var+'AnnuaAvg'] = None

        #         # initiating annualyTotAvg (accumulator variable)
        #         vars(self)[var+'AnnuaTot'] = None

        #         logger.info("Creating the netcdf file for annual average reporting for variable %s.", str(var))

        #         short_name = varDicts.netcdf_short_name[var]
        #         unit       = varDicts.netcdf_unit[var]      
        #         long_name  = varDicts.netcdf_long_name[var]
        #         if long_name == None: long_name = short_name  

        #         # creating netCDF files:
        #         self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                     "_annuaAvg_output.nc",\
        #                                     short_name,unit,long_name)
        # #
        # # -- last day of the year
        # self.outAnnuaEndNC = ["None"]
        # try:
        #     self.outAnnuaEndNC = list(set(self.configuration.reportingOptions['outAnnuaEndNC'].split(",")))
        # except:
        #     pass
        # if self.outAnnuaEndNC[0] != "None":

        #     for var in self.outAnnuaEndNC:

        #         logger.info("Creating the netcdf file for annual end reporting for variable %s.", str(var))

        #         short_name = varDicts.netcdf_short_name[var]
        #         unit       = varDicts.netcdf_unit[var]      
        #         long_name  = varDicts.netcdf_long_name[var]
        #         if long_name == None: long_name = short_name  

        #         # creating netCDF files:
        #         self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                     "_annuaEnd_output.nc",\
        #                                     short_name,unit,long_name)
        # # -- maximum of the year
        # self.outAnnuaMaxNC = ["None"]
        # try:
        #     self.outAnnuaMaxNC = list(set(self.configuration.reportingOptions['outAnnuaMaxNC'].split(",")))
        # except:
        #     pass
        # if self.outAnnuaMaxNC[0] != "None":

        #     for var in self.outAnnuaMaxNC:

        #         logger.info("Creating the netcdf file for annual maximum reporting for variable %s.", str(var))

        #         short_name = varDicts.netcdf_short_name[var]
        #         unit       = varDicts.netcdf_unit[var]      
        #         long_name  = varDicts.netcdf_long_name[var]
        #         if long_name == None: long_name = short_name  

        #         # creating netCDF files:
        #         self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                     "_annuaMax_output.nc",\
        #                                     short_name,unit,long_name)
        
        # list of variables that will be reported:
        self.variables_for_report = self.outDailyTotNC # +\
                                    # self.outMonthTotNC +\
                                    # self.outMonthAvgNC +\
                                    # self.outMonthEndNC +\
                                    # self.outMonthMaxNC +\
                                    # self.outAnnuaTotNC +\
                                    # self.outAnnuaAvgNC +\
                                    # self.outAnnuaEndNC +\
                                    # self.outMonthMaxNC

    def post_processing(self):
        pass                    # see PCR-GLOBWB for details

    def report_forcing_for_debugging(self):
        pass                    # see PCR-GLOBWB for details

    def report_vegetation_phenology_for_debugging(self):
        pass                    # see PCR-GLOBWB for details

    def report_static_maps_for_debugging(self):
        pass                    # see PCR-GLOBWB for details

    def basic_post_processing(self):
        pass                    # see PCR-GLOBWB for details
        
    def additional_post_processing(self):
        pass                    # see PCR-GLOBWB for details
    
    def report(self):

        # # recap all variables
        # self.post_processing()

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
                self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_dailyTot_output.nc",\
                                            short_name,\
                  pcr.pcr2numpy(self.__getattribute__(var),vos.MV),\
                                            timeStamp)

        # # writing monthly output to netcdf files
        # # - cummulative
        # if self.outMonthTotNC[0] != "None":
        #     for var in self.outMonthTotNC:

        #         # introduce variables at the beginning of simulation or
        #         #     reset variables at the beginning of the month
        #         if self._modelTime.timeStepPCR == 1 or \
        #            self._modelTime.day == 1:\
        #            vars(self)[var+'MonthTot'] = pcr.scalar(0.0)

        #         # accumulating
        #         vars(self)[var+'MonthTot'] += vars(self)[var]

        #         # reporting at the end of the month:
        #         if self._modelTime.endMonth == True: 

        #             short_name = varDicts.netcdf_short_name[var]
        #             self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                        "_monthTot_output.nc",\
        #                                        short_name,\
        #               pcr.pcr2numpy(self.__getattribute__(var+'MonthTot'),\
        #                vos.MV),timeStamp)
        # #
        # # - average
        # if self.outMonthAvgNC[0] != "None":
        #     for var in self.outMonthAvgNC:

        #         # only if a accumulator variable has not been defined: 
        #         if var not in self.outMonthTotNC: 

        #             # introduce accumulator at the beginning of simulation or
        #             #     reset accumulator at the beginning of the month
        #             if self._modelTime.timeStepPCR == 1 or \
        #                self._modelTime.day == 1:\
        #                vars(self)[var+'MonthTot'] = pcr.scalar(0.0)

        #             # accumulating
        #             vars(self)[var+'MonthTot'] += vars(self)[var]

        #         # calculating average & reporting at the end of the month:
        #         if self._modelTime.endMonth == True:

        #             vars(self)[var+'MonthAvg'] = vars(self)[var+'MonthTot']/\
        #                                          self._modelTime.day  

        #             short_name = varDicts.netcdf_short_name[var]
        #             self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
        #                                        str(var)+\
        #                                        "_monthAvg_output.nc",\
        #                                        short_name,\
        #               pcr.pcr2numpy(self.__getattribute__(var+'MonthAvg'),\
        #                vos.MV),timeStamp)
        # #
        # # - last day of the month
        # if self.outMonthEndNC[0] != "None":
        #     for var in self.outMonthEndNC:

        #         # reporting at the end of the month:
        #         if self._modelTime.endMonth == True: 

        #             short_name = varDicts.netcdf_short_name[var]
        #             self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
        #                                        str(var)+\
        #                                        "_monthEnd_output.nc",\
        #                                        short_name,\
        #               pcr.pcr2numpy(self.__getattribute__(var),\
        #                vos.MV),timeStamp)
        # #
        # # - maximum
        # if self.outMonthMaxNC[0] != "None":
        #     for var in self.outMonthMaxNC:

        #         # introduce variables at the beginning of simulation or
        #         #     reset variables at the beginning of the month
        #         if self._modelTime.timeStepPCR == 1 or \
        #            self._modelTime.day == 1:\
        #            vars(self)[var+'MonthMax'] = pcr.scalar(0.0)

        #         # find the maximum
        #         vars(self)[var+'MonthMax'] = pcr.max(vars(self)[var], vars(self)[var+'MonthMax'])

        #         # reporting at the end of the month:
        #         if self._modelTime.endMonth == True: 

        #             short_name = varDicts.netcdf_short_name[var]
        #             self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                        "_monthMax_output.nc",\
        #                                        short_name,\
        #               pcr.pcr2numpy(self.__getattribute__(var+'MonthMax'),\
        #                vos.MV),timeStamp)

        # # writing yearly output to netcdf files
        # # - cummulative
        # if self.outAnnuaTotNC[0] != "None":
        #     for var in self.outAnnuaTotNC:

        #         # introduce variables at the beginning of simulation or
        #         #     reset variables at the beginning of the year
        #         if self._modelTime.timeStepPCR == 1 or \
        #            self._modelTime.doy == 1:\
        #            vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)

        #         # accumulating
        #         vars(self)[var+'AnnuaTot'] += vars(self)[var]

        #         # reporting at the end of the year:
        #         if self._modelTime.endYear == True: 

        #             short_name = varDicts.netcdf_short_name[var]
        #             self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
        #                                        str(var)+\
        #                                        "_annuaTot_output.nc",\
        #                                        short_name,\
        #               pcr.pcr2numpy(self.__getattribute__(var+'AnnuaTot'),\
        #                vos.MV),timeStamp)

        # # - average
        # if self.outAnnuaAvgNC[0] != "None":
        #     for var in self.outAnnuaAvgNC:

        #         # only if a accumulator variable has not been defined: 
        #         if var not in self.outAnnuaTotNC: 

        #             # introduce accumulator at the beginning of simulation or
        #             #     reset accumulator at the beginning of the year
        #             if self._modelTime.timeStepPCR == 1 or \
        #                self._modelTime.doy == 1:\
        #                vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)

        #             # accumulating
        #             vars(self)[var+'AnnuaTot'] += vars(self)[var]

        #         # calculating average & reporting at the end of the year:
        #         if self._modelTime.endYear == True:

        #             vars(self)[var+'AnnuaAvg'] = vars(self)[var+'AnnuaTot']/\
        #                                          self._modelTime.doy  

        #             short_name = varDicts.netcdf_short_name[var]
        #             self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
        #                                        str(var)+\
        #                                        "_annuaAvg_output.nc",\
        #                                        short_name,\
        #               pcr.pcr2numpy(self.__getattribute__(var+'AnnuaAvg'),\
        #                vos.MV),timeStamp)
        # #
        # # -last day of the year
        # if self.outAnnuaEndNC[0] != "None":
        #     for var in self.outAnnuaEndNC:

        #         # calculating average & reporting at the end of the year:
        #         if self._modelTime.endYear == True:

        #             short_name = varDicts.netcdf_short_name[var]
        #             self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
        #                                        str(var)+\
        #                                        "_annuaEnd_output.nc",\
        #                                        short_name,\
        #               pcr.pcr2numpy(self.__getattribute__(var),\
        #                vos.MV),timeStamp)
        # #
        # # - maximum
        # if self.outAnnuaMaxNC[0] != "None":
        #     for var in self.outAnnuaMaxNC:

        #         # introduce variables at the beginning of simulation or
        #         #     reset variables at the beginning of the year
        #         if self._modelTime.timeStepPCR == 1 or \
        #            self._modelTime.doy == 1:\
        #            vars(self)[var+'AnnuaMax'] = pcr.scalar(0.0)

        #         # find the maximum
        #         vars(self)[var+'AnnuaMax'] = pcr.max(vars(self)[var], vars(self)[var+'AnnuaMax'])

        #         # reporting at the end of the year:
        #         if self._modelTime.endYear == True: 

        #             short_name = varDicts.netcdf_short_name[var]
        #             self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
        #                                     str(var)+\
        #                                        "_annuaMax_output.nc",\
        #                                        short_name,\
        #               pcr.pcr2numpy(self.__getattribute__(var+'AnnuaMax'),\
        #                vos.MV),timeStamp)
       
