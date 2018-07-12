#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
import sys

from pcraster.framework import DynamicModel
from pcraster.framework import DynamicFramework

from Configuration import Configuration
from CurrTimeStep import ModelTime
from Reporting import *
# from spinUp import SpinUp

from AquaCrop import AquaCrop
from FAO56 import FAO56

import logging
logger = logging.getLogger(__name__)

class DeterministicRunner(DynamicModel):

    # def __init__(self, configuration, modelTime, initialState = None):
    #     DynamicModel.__init__(self)
    #     self.modelTime = modelTime
    #     self.model = AquaCrop(configuration, modelTime, initialState)
    #     self.reporting = Reporting_AquaCrop(configuration, self.model, modelTime)

    def initial(self):
        pass

    def dynamic(self):

        # recalculate current model time
        self.modelTime.update(self.currentTimeStep())

        # update model
        self.model.dynamic()

        # do any reporting
        self.reporting.report()

class run_AquaCrop(DeterministicRunner):
    def __init__(self, configuration, modelTime, initialState = None):
        DynamicModel.__init__(self)
        self.modelTime = modelTime
        self.model = AquaCrop(configuration, modelTime, initialState)
        self.model.initial()
        self.reporting = Reporting(configuration, self.model, modelTime)

class run_FAO56(DeterministicRunner):
    def __init__(self, configuration, modelTime, initialState = None):
        DynamicModel.__init__(self)
        self.modelTime = modelTime
        self.model = FAO56(configuration, modelTime, initialState)
        self.model.initial()
        self.reporting = Reporting(configuration, self.model, modelTime)

def main():

    # TODO: print disclaimer
    # disclaimer.print_disclaimer()

    # get the full path of the configuration/ini file provided
    # as system argument
    iniFileName = os.path.abspath(sys.argv[1])
    model = sys.argv[2].lower()

    # TODO: debug option
    debug_mode = False
    # if len(sys.argv) > 2:
    #     if (sys.argv[2] == "debug": debug_mode = True

    # object to handle configuration/ini file
    configuration = Configuration(iniFileName=iniFileName, debug_mode=debug_mode)

    # timestep information
    currTimeStep = ModelTime()

    # NB spin-up currently not implemented (nor is it
    # implemented in AquaCrop-OS)
    initial_state = None

    currTimeStep.getStartEndTimeSteps(configuration.globalOptions['startTime'], configuration.globalOptions['endTime'])
    currTimeStep.update(1)      # this essentially allows us to call read_forcings in AquaCrop.__init__() method
    
    logger.info('Transient simulation run has started')
    deterministic_runner = None
    if (model == 'aquacrop'):
        deterministic_runner = run_AquaCrop(configuration, currTimeStep, initial_state)
    elif model == 'fao56':
        deterministic_runner = run_FAO56(configuration, currTimeStep, initial_state)
    # TODO: error handling
    
    dynamic_framework = DynamicFramework(deterministic_runner, currTimeStep.nrOfTimeSteps)
    dynamic_framework.setQuiet(True)
    dynamic_framework.run()

if __name__ == '__main__':
    # disclaimer.print_disclaimer(with_logger = True)
    sys.exit(main())
