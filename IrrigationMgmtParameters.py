#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# AquaCrop
#
import os
import numpy as np
import pcraster as pcr
import VirtualOS as vos
import netCDF4 as nc
import datetime as datetime
import calendar as calendar

class IrrigationMgmtParameters(object):

    def __init__(self, IrrigationMgmtParameters_variable):
        self.var = IrrigationMgmtParameters_variable

    def initial(self):
        self.var.irrMgmtParameterFileNC = self.var._configuration.irrMgmtOptions['irrMgmtParameterNC']
        self.var.parameter_names = ['IrrMethod','IrrInterval','SMT1','SMT2','SMT3','SMT4','MaxIrr','AppEff','NetIrrSMT','WetSurf']
        for var in self.var.parameter_names:
            nm = '_' + var
            vars(self.var)[nm] = vos.netcdf2PCRobjCloneWithoutTime(
                self.var.irrMgmtParameterFileNC,
                var,
                cloneMapFileName=self.var.cloneMap)

        # check if an irrigation schedule file is required
        if np.sum(self.var._IrrMethod == 3) > 0:
            if self.var._configuration.irrMgmtOptions['irrScheduleNC'] != "None":
                self.var._configuration.irrMgmtOptions['irrScheduleNC'] = vos.getFullPath(self.var._configuration.irrMgmtOptions[item], self.var._configuration.globalOptions['inputDir'])
                self.var.irrScheduleFileNC = self.var._configuration.irrMgmtOptions['irrScheduleNC']
            else:
                logger.error('IrrMethod equals 3 in some or all places, but irrScheduleNC is not set in configuration file')

        else:
            self.var.irrScheduleFileNC = None

    def dynamic(self):
        I,J,K = np.ogrid[:self.var.nRotation,:self.var.nLat,:self.var.nLon]
        for var in self.var.parameter_names:
            vars(self.var)[var] = getattr(self.var, '_' + var)[self.var.CropIndex,I,J,K]
