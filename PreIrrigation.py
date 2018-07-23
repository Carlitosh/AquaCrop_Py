#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class PreIrrigation(object):
    """Class to represent pre-irrigation when in net irrigation 
    mode.
    """
    def __init__(self, PreIrrigation_variable):
        self.var = PreIrrigation_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.PreIrr = np.copy(arr_zeros)
        self.var.IrrNet = np.copy(arr_zeros)
        
    def dynamic(self):
        # Expand dz and dzsum to crop, lat, lon
        arr_ones = np.ones((self.var.nCrop, self.var.nLat, self.var.nLon))[:,None,:,:]
        dz = (self.var.dz[None,:,None,None] * arr_ones)
        dzsum = (self.var.dzsum[None,:,None,None] * arr_ones)

        # Calculate pre-irrigation requirement
        rootdepth = np.maximum(self.var.Zmin, self.var.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        netirrsmt = self.var.NetIrrSMT[:,None,:,:] * np.ones((self.var.nComp))[None,:,None,None]
        thCrit = (self.var.th_wp_comp + ((netirrsmt / 100) * (self.var.th_fc_comp - self.var.th_wp_comp)))

        # Conditions for applying pre-irrigation
        cond0 = ((self.var.IrrMethod == 4) & (self.var.DAP == 1))
        cond1 = (np.broadcast_to(cond0[:,None,:,:], self.var.th.shape)
                 & ((dzsum - dz) < np.broadcast_to(rootdepth[:,None,:,:], self.var.th.shape))
                 & (self.var.th < thCrit))

        # Update pre-irrigation and root zone water content (mm)
        PreIrr_req = ((thCrit - self.var.th) * 1000 * dz)
        PreIrr_req[np.logical_not(cond1)] = 0
        self.var.PreIrr = np.sum(PreIrr_req, axis=1)
                
    def add_pre_irrigation(self):
        self.var.IrrNet += self.var.PreIrr
