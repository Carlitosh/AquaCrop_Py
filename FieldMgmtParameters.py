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

    # TODO: make sure field management parameters have dimensions nRotation, nLat, nLon
    def __init__(self, FieldMgmtParameters_variable):
        self.var = FieldMgmtParameters_variable
        
    def initial(self):
        self.var.fieldMgmtParameterFileNC = self.var._configuration.fieldMgmtOptions['fieldMgmtParameterNC']
        self.var.parameter_names = ['Mulches','MulchPctGS','MulchPctOS','fMulch','Bunds','zBund','BundWater']
        for var in self.var.parameter_names:
            # nm = '_' + var
            vars(self.var)[var] = vos.netcdf2PCRobjCloneWithoutTime(
                self.var.fieldMgmtParameterFileNC,
                var,
                cloneMapFileName=self.var.cloneMap)

    def dynamic(self):
        pass
