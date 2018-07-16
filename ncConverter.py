#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
import sys
import datetime
import time
import re
import glob
import subprocess
import netCDF4 as nc
import numpy as np
import pcraster as pcr
import VirtualOS as vos

# TODO: defined the dictionary (e.g. filecache = dict()) to avoid open and closing files

class np2netCDF():
    
    def __init__(self,configuration,model,specificAttributeDictionary=None):

	# Set clone map
        pcr.setclone(configuration.cloneMap)
        cloneMap = pcr.boolean(1.0)  # map with all cell values equal to 1

        # Retrieve latitudes and longitudes from clone map
        self.latitudes  = np.unique(pcr.pcr2numpy(pcr.ycoordinate(cloneMap), vos.MV))[::-1]
        self.longitudes = np.unique(pcr.pcr2numpy(pcr.xcoordinate(cloneMap), vos.MV))
        self.crops  = np.arange(1, model.nCrop + 1)
        self.depths = np.arange(1, model.nComp + 1)
        
        # Let users decide what their preference regarding latitude order
        self.netcdf_y_orientation_follow_cf_convention = False
        if 'netcdf_y_orientation_follow_cf_convention' in configuration.reportingOptions.keys() and\
            configuration.reportingOptions['netcdf_y_orientation_follow_cf_convention'] == "True":
            msg = "Latitude (y) orientation for output netcdf files start from the bottom to top."
            self.netcdf_y_orientation_follow_cf_convention = True
            self.latitudes  = np.unique(pcr.pcr2numpy(pcr.ycoordinate(cloneMap), vos.MV))
        
        # Set general netcdf attributes (based on the information given in the ini/configuration file) 
        self.set_general_netcdf_attributes(configuration, specificAttributeDictionary)
        
        # netcdf format and zlib setup 
        self.format = 'NETCDF3_CLASSIC'
        self.zlib = False
        if "formatNetCDF" in configuration.reportingOptions.keys():
            self.format = str(configuration.reportingOptions['formatNetCDF'])
        if "zlib" in configuration.reportingOptions.keys():
            if configuration.reportingOptions['zlib'] == "True": self.zlib = True

        # # if given in the ini file, use the netcdf as given in the section 'specific_attributes_for_netcdf_output_files'
        # if 'specific_attributes_for_netcdf_output_files' in configuration.allSections:
        #     for key in configuration.specific_attributes_for_netcdf_output_files.keys():

        #         self.attributeDictionary[key] = configuration.specific_attributes_for_netcdf_output_files[key]
                
        #         if self.attributeDictionary[key] == "None": self.attributeDictionary[key] = ""

        #         if key == "history" and self.attributeDictionary[key] == "Default":
        #             self.attributeDictionary[key] = \
        #                             'created on ' + datetime.datetime.today().isoformat(' ')
        #         if self.attributeDictionary[key] == "Default" and\
        #           (key == "date_created" or key == "date_issued"):
        #             self.attributeDictionary[key] = datetime.datetime.today().isoformat(' ')
                    
    def set_general_netcdf_attributes(self,configuration,specificAttributeDictionary=None):
        """Function to set general netCDF attributes"""
        
        # netCDF attributes (based on the configuration file or specificAttributeDictionary):
        self.attributeDictionary = {}
        if specificAttributeDictionary == None:
            self.attributeDictionary['institution'] = configuration.globalOptions['institution']
            self.attributeDictionary['title'      ] = configuration.globalOptions['title'      ]
            self.attributeDictionary['description'] = configuration.globalOptions['description']
        else:
            self.attributeDictionary['institution'] = specificAttributeDictionary['institution']
            self.attributeDictionary['title'      ] = specificAttributeDictionary['title'      ]
            self.attributeDictionary['description'] = specificAttributeDictionary['description']

    def createNetCDF(self, ncFileName, varName, varUnits, varDims, longName = None):
        """Function to create output netCDF"""
        
        rootgrp = nc.Dataset(ncFileName,'w',format= self.format)

        # Create dimensions - time is unlimited, others are fixed
        if 'crop' in varDims: rootgrp.createDimension('crop',len(self.crops))
        if 'time' in varDims:     rootgrp.createDimension('time',None)
        if 'depth' in varDims:    rootgrp.createDimension('depth', len(self.depths))
        if 'lat' in varDims:      rootgrp.createDimension('lat',len(self.latitudes))
        if 'lon' in varDims:      rootgrp.createDimension('lon',len(self.longitudes))

        # define variables (i4 = 32-bit integer, f4 = 32-bit floating point)
        if 'crop' in varDims:
            crop = rootgrp.createVariable('crop','i4',('crop',))
            crop.standard_name = 'crop'
            crop.long_name = 'crop'
            crop[:] = self.crops

        if 'time' in varDims:
            date_time = rootgrp.createVariable('time','f4',('time',))
            date_time.standard_name = 'time'
            date_time.long_name = 'Days since 1901-01-01'
            date_time.units = 'Days since 1901-01-01' 
            date_time.calendar = 'standard'

        # if includeDepthDimension:
        if 'depth' in varDims:
            depth = rootgrp.createVariable('depth','f4',('depth',))
            depth.standard_name = 'depth'
            depth.long_name = 'depth'
            depth.units = 'meter'
            depth.positive = 'down'
            depth[:] = self.depths

        if 'lat' in varDims:
            lat = rootgrp.createVariable('lat','f4',('lat',))
            lat.long_name = 'latitude'
            lat.units = 'degrees_north'
            lat.standard_name = 'latitude'
            lat[:]= self.latitudes

        if 'lon' in varDims:
            lon = rootgrp.createVariable('lon','f4',('lon',))
            lon.long_name = 'longitude'
            lon.units = 'degrees_east'
            lon.standard_name = 'longitude'
            lon[:]= self.longitudes

        dims = varDims

        # Add variable to NetCDF
        shortVarName = varName
        longVarName  = varName
        if longName != None: longVarName = longName
        var = rootgrp.createVariable(shortVarName,
                                     'f4',
                                     dims,
                                     fill_value=vos.MV,
                                     zlib=self.zlib)
        var.standard_name = varName
        var.long_name = longVarName
        var.units = varUnits

        attributeDictionary = self.attributeDictionary
        for k, v in attributeDictionary.items(): setattr(rootgrp,k,v)
        rootgrp.sync()
        rootgrp.close()
                
    def data2NetCDF(self, ncFileName, shortVarName, dims, varField, timeStamp, posCnt = None):
        """Function to write data to netCDF"""

        rootgrp = nc.Dataset(ncFileName,'a')

        date_time = rootgrp.variables['time']
        if posCnt == None: posCnt = len(date_time)
        date_time[posCnt] = nc.date2num(timeStamp,date_time.units,date_time.calendar)

        # flip variable if necessary (to follow cf_convention)
        if self.netcdf_y_orientation_follow_cf_convention: varField = np.flipud(varField)

        if 'depth' in dims:
            rootgrp.variables[shortVarName][:,posCnt,:,:,:] = varField
        else:
            if 'crop' in dims:
                rootgrp.variables[shortVarName][:,posCnt,:,:] = varField
            else:
                rootgrp.variables[shortVarName][posCnt,:,:] = varField
                        
        rootgrp.sync()
        rootgrp.close()

    def close(self, ncFileName):
        """Function to close netCDF file"""
        rootgrp = nc.Dataset(ncFileName,'w')
        rootgrp.close()
