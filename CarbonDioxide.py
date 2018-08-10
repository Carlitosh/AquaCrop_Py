#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import VirtualOS as vos

import logging
logger = logging.getLogger(__name__)

class CarbonDioxide(object):

    def __init__(self, CarbonDioxide_variable):
        self.var = CarbonDioxide_variable

    def initial(self):
        # reference concentration
        self.var.RefConc = 369.41        
        self.var.co2FileNC = self.var._configuration.carbonDioxideOptions['carbonDioxideNC']

        # variable names      
        self.var.co2VarName = 'co2' 
        if 'co2VariableName' in self.var._configuration.carbonDioxideOptions:
            self.var.co2VarName = self.var._configuration.carbonDioxideOptions['co2VariableName']
        self.var.co2_set_per_year  = False
        
    def dynamic(self):

        # TODO: only read data if the year has changed
        
        if self.var._modelTime.timeStepPCR == 1 or self.var._modelTime.doy == 1:
            date = '%04i-%02i-%02i' %(self.var._modelTime.year, 1, 1)
            self.var.conc = vos.netcdf2PCRobjClone(self.var.co2FileNC,
                                               self.var.co2VarName,
                                               date,
                                               useDoy = None,
                                               cloneMapFileName = self.var.cloneMap,
                                               LatitudeLongitude = True)
