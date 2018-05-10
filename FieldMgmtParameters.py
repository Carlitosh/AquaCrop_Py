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

class FieldMgmtParameters(object):

    def __init__(self, iniItems, landmask):
        object.__init__(self)

        self.cloneMap = iniItems.cloneMap
        self.landmask = landmask

        attr = vos.getMapAttributesALL(self.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])

        self.fieldMgmtParameterFileNC = iniItems.fieldMgmtOptions['fieldMgmtParameterNC']

    def read(self):
        """Function to read field management input parameters"""
        
        self.parameter_names = ['Mulches','MulchPctGS','MulchPctOS','fMulch','Bunds','zBund','BundWater']
        for var in self.parameter_names:
            vars(self)[var] = vos.netcdf2PCRobjCloneWithoutTime(
                self.fieldMgmtParameterFileNC,
                var,
                cloneMapFileName=self.cloneMap)            
