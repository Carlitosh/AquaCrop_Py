#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

# The purpose of this file is to parse the configuration file

import ConfigParser
import os
import sys
import virtualOS as vos
import time
import datetime
import shutil
import glob

import logging
logger = logging.getLogger(__name__)

class Configuration(object):

    def __init__(self, iniFileName, debug_mode = False, no_modification = True, system_arguments = None, relative_ini_meteo_paths = False):
        object.__init__(self)

        if iniFileName is None:
            raise Exception('Error: No configuration file specified')
        self._timestamp = datetime.datetime.now()

        print self._timestamp

        # get the full path of iniFileName
        self.iniFileName = os.path.abspath(iniFileName)

        # debug option
        self.debug_mode = debug_mode

        # save cwd
        self._cwd = os.getcwd()

        # read configuration file
        self.parse_configuration_file(self.iniFileName)

        # set configuration
        if no_modification: self.set_configuration(system_arguments)

        # set the main output directory
        self.main_output_directory = self.globalOptions['outputDir']

    def set_configuration(self, system_arguments = None):

        # set all paths
        self.set_input_files()
        self.create_output_directories()

        # initialize logging 
        self.initialize_logging("Default", system_arguments)

        # copy ini file
        self.backup_configuration()

        # repair key names of initial conditions
        self.repair_ini_key_names()

    def initialize_logging(self, log_file_location = "Default", system_arguments = None):
        """
        Initialize logging. Prints to both the console and a log file, at configurable levels
        """

        # set root logger to debug level        
        logging.getLogger().setLevel(logging.DEBUG)

        # logging format 
        formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')

        # default logging levels
        log_level_console    = "INFO"
        log_level_file       = "INFO"
        # order: DEBUG, INFO, WARNING, ERROR, CRITICAL
        
        # log level based on ini/configuration file:
        if "log_level_console" in self.globalOptions.keys():
            log_level_console = self.globalOptions['log_level_console']        
        if "log_level_file" in self.globalOptions.keys():
            log_level_file = self.globalOptions['log_level_file']        

        # log level for debug mode:
        if self.debug_mode == True: 
            log_level_console = "DEBUG"
            log_level_file    = "DEBUG"

        console_level = getattr(logging, log_level_console.upper(), logging.INFO)
        if not isinstance(console_level, int):
            raise ValueError('Invalid log level: %s', log_level_console)
        
        # create handler, add to root logger
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        console_handler.setLevel(console_level)
        logging.getLogger().addHandler(console_handler)

        # log file name (and location)
        if log_file_location != "Default":  self.logFileDir = log_file_location
        log_filename = self.logFileDir + os.path.basename(self.iniFileName) + '_' + str(self._timestamp.isoformat()).replace(":",".") + '.log'

        file_level = getattr(logging, log_level_file.upper(), logging.DEBUG)
        if not isinstance(console_level, int):
            raise ValueError('Invalid log level: %s', log_level_file)

        # create handler, add to root logger
        file_handler = logging.FileHandler(log_filename)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(file_level)
        logging.getLogger().addHandler(file_handler)
        
        # file name for debug log 
        dbg_filename = self.logFileDir + os.path.basename(self.iniFileName) + '_' +  str(self._timestamp.isoformat()).replace(":",".") + '.dbg'

        # create handler, add to root logger
        debug_handler = logging.FileHandler(dbg_filename)
        debug_handler.setFormatter(formatter)
        debug_handler.setLevel(logging.DEBUG)
        logging.getLogger().addHandler(debug_handler)

        # # print disclaimer
        # disclaimer.print_disclaimer(with_logger = True)
        
        logger.info('Model run started at %s', self._timestamp)
        logger.info('Logging output to %s', log_filename)
        logger.info('Debugging output to %s', dbg_filename)
        
        if system_arguments != None:
            logger.info('The system arguments given to execute this run: %s', system_arguments)
        

    def backup_configuration(self):

        # copy ini file to logDir
        shutil.copy(self.iniFileName, self.logFileDir + os.path.basename(self.iniFileName) + '_' + str(self._timestamp.isoformat()).replace(":",".") + '.ini')
        
        
    def parse_configuration_file(self, modelFileName):
        config = ConfigParser.ConfigParser()
        config.optionxform = str
        config.read(modelFileName)

        # sections in configuration file
        self.allSections = config.sections()

        # read all sections
        for section in self.allSections:
            vars(self)[section] = {}
            options = config.options(section)
            for option in options:
                val = config.get(section, option)
                self.__getattribute__(section)[option] = val

    def set_input_files(self):
        
        # set clone map
        self.cloneMap = vos.getFullPath(self.globalOptions['cloneMap'], self.globalOptions['inputDir'])

        # meteorological input files
        meteoInputFiles = ['precipitationNC','temperatureNC','refETPotFileNC']
        for item in meteoInputFiles:
            if self.meteoOptions[item] != "None":
                self.meteoOptions[item] = vos.getFullPath(self.meteoOptions[item], self.globalOptions['inputDir'])

        # CO2 concentration input file
        co2InputFiles = ['carbonDioxideNC']
        for item in co2InputFiles:
            if self.carbonDioxideOptions[item] != "None":
                self.carbonDioxideOptions[item] = vos.getFullPath(self.carbondDioxideOptions[item], self.globalOptions['inputDir'])
                
        # groundwater input file
        item = 'groundwaterNC'
        if self.groundwaterOptions[item] != "None":
            self.groundwaterOptions[item] = vos.getFullPath(self.groundwaterOptions[item], self.globalOptions['inputDir'])

        # soil parameter input file
        item = 'soilAndTopoNC'
        if self.soilOptions[item] != "None":
            self.soilOptions[item] = vos.getFullPath(self.soilOptions[item], self.globalOptions['inputDir'])

        # crop parameter input file
        item = 'cropParameterNC'
        if self.cropOptions[item] != "None":
            self.cropOptions[item] = vos.getFullPath(self.cropOptions[item], self.globalOptions['inputDir'])

        # rotation management file
        item = 'rotationParameterNC'
        if self.rotationOptions[item] != "None":
            self.rotationOptions[item] = vos.getFullPath(self.rotationOptions[item], self.globalOptions['inputDir'])
            
        # irrigation management input file
        irrMgmtInputFiles = ['irrMgmtParameterNC']
        for item in irrMgmtInputFiles:
            if self.irrMgmtOptions[item] != "None":
                self.irrMgmtOptions[item] = vos.getFullPath(self.irrMgmtOptions[item], self.globalOptions['inputDir'])
            
        # field management input file
        item = 'fieldMgmtParameterNC'
        if self.fieldMgmtOptions[item] != "None":
            self.fieldMgmtOptions[item] = vos.getFullPath(self.fieldMgmtOptions[item], self.globalOptions['inputDir'])
            
        # land surface initial conditions file
        item = 'initialConditionNC'
        if self.landSurfaceOptions[item] != "None":
            self.landSurfaceOptions[item] = vos.getFullPath(self.landSurfaceOptions[item], self.globalOptions['inputDir'])
            
    def create_output_directories(self):

        cleanOutputDir = False
        if cleanOutputDir:
            try:
                shutil.rmtree(self.globalOptions['outputDir'])
            except:
                pass

        try:
            os.makedirs(self.globalOptions['outputDir'])
        except:
            pass

        # make temp directory
        self.tmpDir = vos.getFullPath("tmp/", self.globalOptions['outputDir'])

        if os.path.exists(self.tmpDir):
            shutil.rmtree(self.tmpDir)

        os.makedirs(self.tmpDir)

        # make netcdf directory
        self.outNCDir = vos.getFullPath("netcdf/", self.globalOptions['outputDir'])
        
        if os.path.exists(self.outNCDir):
            shutil.rmtree(self.outNCDir)

        os.makedirs(self.outNCDir)

        # make and populate backup directory for Python
        # scripts

        # SM: I guess the idea here is that the output
        # directory contains the exact version of the model
        # used to create the output
        self.scriptDir = vos.getFullPath("scripts/", self.globalOptions['outputDir'])

        if os.path.exists(self.scriptDir):
            shutil.rmtree(self.scriptDir)

        os.makedirs(self.scriptDir)

        # working/starting directory where all Python scripts
        # are located
        path_of_this_module = os.path.abspath(os.path.dirname(__file__))
        self.starting_directory = path_of_this_module

        all_files = glob.glob(os.path.join(path_of_this_module, '*.py'))
        for filename in all_files:
            shutil.copy(filename, self.scriptDir)

        # make log directory
        self.logFileDir = vos.getFullPath("log/", self.globalOptions['outputDir'])

        cleanLogDir = True
        if os.path.exists(self.logFileDir) and cleanLogDir:
            shutil.rmtree(self.logFileDir)

        os.makedirs(self.logFileDir)

        # make end state directory
        self.endStateDir = vos.getFullPath("states/", self.globalOptions['outputDir'])

        if os.path.exists(self.endStateDir):
            shutil.rmtree(self.endStateDir)

        os.makedirs(self.endStateDir)

        # PCRaster maps directory? I don't think this is needed
        
    def repair_ini_key_names(self):
        """
        If needed, change/modify some key names for initial condition fields. 
        This is introduced because Edwin was very stupid as once he changed some key names of initial conditions! Yet, it is also useful particularly for runs without complete ini files. 
        """

        # temporal resolution of the model
        self.timeStep = 1.0
        self.timeStepUnit = "day"
        if 'timeStep' in self.globalOptions.keys() and \
           'timeStepUnit'in self.globalOptions.keys():

            if float(self.globalOptions['timeStep']) != 1.0 or \
                     self.globalOptions['timeStepUnit'] != "day":
                logger.error('The model runs only on daily time step. Please check your ini/configuration file')
                self.timeStep     = None
                self.timeStepUnit = None

        # TODO: additional ini key names
        
        # groundwater options
        # ===================
        
        if 'WaterTable' not in self.groundwaterOptions.keys():
            # TODO: add warning to say that groundwater is not considered
            # TODO: add check that groundwater file is supplied if WaterTable = 1
            self.groundwaterOptions['WaterTable'] = "0"
        
        if self.groundwaterOptions['WaterTable'] == "0":
            watertable = False
        elif self.groundwaterOptions['WaterTable'] == "1":
            watertable = True
        else:
            # TODO: print warning ("e.g. WaterTable must equal 1 or 0")
            self.groundwaterOptions['WaterTable'] == "0"
            watertable = False
                 

        if watertable and 'VariableWaterTable' not in self.groundwaterOptions.keys():
            # TODO: add warning (also, is this appropriate?)
            self.groundwaterOptions['WaterTable'] = "0"
            watertable = False

        if watertable and 'groundwaterVariableName' not in self.groundwaterOptions.keys():
            # TODO: add warning (also, is this appropriate?)
            self.groundwaterOptions['WaterTable'] = "0"
            watertable = False

        if watertable and 'groundwaterNC' not in self.groundwaterOptions.keys():
            # TODO: add warning
            self.groundwaterOptions['WaterTable'] = "0"
            watertable = False
            
        # irrigation options
        # ==================
        
        # irrigation schedule file 
        if 'irrScheduleNC' not in self.irrMgmtOptions.keys():
            self.irrMgmtOptions['irrScheduleNC'] = "None"
