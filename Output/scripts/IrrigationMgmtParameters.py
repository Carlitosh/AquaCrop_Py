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

    def __init__(self, iniItems, landmask):
        object.__init__(self)

        self.cloneMap = iniItems.cloneMap
        self.landmask = landmask

        attr = vos.getMapAttributesALL(self.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])

        self.irrMgmtParameterFileNC = iniItems.irrMgmtOptions['irrMgmtParameterNC']

    def read(self):
        """Function to read irrigation management input parameters"""
        
        self.parameter_names = ['IrrMethod','IrrInterval','SMT1','SMT2','SMT3','SMT4','MaxIrr','AppEff','NetIrrSMT','WetSurf']
        for var in self.parameter_names:
            vars(self)[var] = vos.netcdf2PCRobjCloneWithoutTime(
                self.irrMgmtParameterFileNC,
                var,
                cloneMapFileName=self.cloneMap)

        # # create array of soil moisture target for the respective growth stages
        # self.SMT = np.concatenate((self.SMT1[None,:], self.SMT2[None,:], self.SMT3[None,:], self.SMT4[None,:]), axis=0)

        # check if an irrigation schedule file is required
        if np.sum(self.IrrMethod == 3) > 0:
            if self.irrMgmtOptions['irrScheduleNC'] != "None":
                self.irrMgmtOptions['irrScheduleNC'] = vos.getFullPath(self.irrMgmtOptions[item], self.globalOptions['inputDir'])
                self.irrScheduleFileNC = irrMgmtOptions['irrScheduleNC']
            else:
                logger.error('IrrMethod equals 3 in some or all places, but irrScheduleNC is not set in configuration file')

        else:
            self.irrScheduleFileNC = None
