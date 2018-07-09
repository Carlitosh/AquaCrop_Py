#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from hydrology_funs import *

import logging
logger = logging.getLogger(__name__)

class Infiltration(object):
    """Class to infiltrate incoming water"""
    
    def __init__(self, Infiltration_variable):
        self.var = Infiltration_variable

    def initial(self):        
        pass
    
    def dynamic(self):
        """Function to infiltrate incoming water (rainfall and 
        irrigation)
        """
        # Update infiltration rate for irrigation
        self.var.Infl += self.var.Irr * (self.var.AppEff / 100)

        SurfaceStorage, RunoffIni, ToStore = update_surface_storage(
            self.var.Bunds,
            self.var.zBund,
            self.var.ksat_comp,
            self.var.Infl,
            self.var.SurfaceStorage)
        
        Runoff, FluxOut, DeepPerc, thnew = infiltration(
            ToStore,
            self.var.FluxOut,
            self.var.th,
            self.var.th_s_comp,
            self.var.th_fc_comp,
            self.var.th_fc_adj,
            self.var.tau_comp,
            self.var.ksat_comp,
            self.var.dz)
        
        # Update total runoff
        Runoff += RunoffIni

        # Update surface storage (if bunds are present)
        self.var.SurfaceStorage = SurfaceStorage
        cond5 = ((Runoff > RunoffIni) & (self.var.Bunds == 1) & self.var.zBund > 0.001)
        self.var.SurfaceStorage[cond5] += (Runoff - RunoffIni)[cond5]

        # Limit surface storage to bund height: additional water above top of
        # bunds becomes runoff, and surface storage equals bund height
        cond51 = (cond5 & (self.var.SurfaceStorage > (self.var.zBund * 1000)))
        Runoff[cond51] = (RunoffIni + (self.var.SurfaceStorage - (self.var.zBund * 1000)))[cond51]
        self.var.SurfaceStorage[cond51] = (self.var.zBund * 1000)[cond51]
        cond52 = (cond5 & np.logical_not(cond51))
        Runoff[cond52] = RunoffIni[cond52]

        # Update water content, deep percolation, surface runoff, infiltration
        self.var.th = np.copy(thnew)
        self.var.DeepPerc += DeepPerc
        self.var.Infl -= Runoff
        self.var.Runoff += Runoff
        
        # print self.var.th[0,0,0,0]
