#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class RootZoneWater(object):
    def __init__(self, RootZoneWater_variable):
        self.var = RootZoneWater_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLon, self.var.nLat))
        self.thRZ_Act  = np.copy(arr_zeros)
        self.thRZ_Sat  = np.copy(arr_zeros)
        self.thRZ_Fc   = np.copy(arr_zeros)
        self.thRZ_Wp   = np.copy(arr_zeros)
        self.thRZ_Dry  = np.copy(arr_zeros)
        self.thRZ_Aer  = np.copy(arr_zeros)
        self.TAW       = np.copy(arr_zeros)
        self.Dr        = np.copy(arr_zeros)
        self.Wr        = np.copy(arr_zeros)
        
    def dynamic(self):
        """Function to calculate actual and total available water in the 
        root zone at current time step
        """
        dims = self.var.th.shape
        arr_ones = np.ones((dims[1], dims[2], dims[3]))
        dz = self.var.dz[:,None,None,None] * arr_ones
        dzsum = self.var.dzsum[:,None,None,None] * arr_ones

        # Calculate root zone water content and available water
        rootdepth = np.maximum(self.var.Zmin, self.var.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        comp_sto = (np.round((dzsum - dz) * 1000) < np.round(rootdepth * 1000))

        # Fraction of compartment covered by root zone (zero in compartments
        # NOT covered by the root zone)
        factor = 1 - ((dzsum - rootdepth) / dz)
        factor = np.clip(factor, 0, 1)
        factor[np.logical_not(comp_sto)] = 0

        # Water storages in root zone (mm) - initially compute value in each
        # compartment, then sum to get overall root zone storages
        Wr_comp = factor * 1000 * self.var.th * dz
        WrS_comp = factor * 1000 * self.var.th_s_comp * dz
        WrFC_comp = factor * 1000 * self.var.th_fc_comp * dz
        WrWP_comp = factor * 1000 * self.var.th_wp_comp * dz
        WrDry_comp = factor * 1000 * self.var.th_dry_comp * dz

        # Water storage in root zone at aeration stress threshold (mm)
        WrAer_comp = factor * 1000 * (self.var.th_s_comp - (self.var.Aer / 100)) * dz

        Wr = np.sum(Wr_comp, axis=0)
        Wr[Wr < 0] = 0
        WrS = np.sum(WrS_comp, axis=0)
        WrFC = np.sum(WrFC_comp, axis=0)
        WrWP = np.sum(WrWP_comp, axis=0)
        WrDry = np.sum(WrDry_comp, axis=0)
        WrAer = np.sum(WrAer_comp, axis=0)

        # Convert depths to m3/m3
        self.var.thRZ_Act = np.divide(Wr, rootdepth * 1000, out=np.zeros_like(Wr), where=rootdepth!=0)
        self.var.thRZ_Sat = np.divide(WrS, rootdepth * 1000, out=np.zeros_like(WrS), where=rootdepth!=0)
        self.var.thRZ_Fc  = np.divide(WrFC, rootdepth * 1000, out=np.zeros_like(WrFC), where=rootdepth!=0)
        self.var.thRZ_Wp  = np.divide(WrWP, rootdepth * 1000, out=np.zeros_like(WrWP), where=rootdepth!=0)
        self.var.thRZ_Dry = np.divide(WrDry, rootdepth * 1000, out=np.zeros_like(WrDry), where=rootdepth!=0)
        self.var.thRZ_Aer = np.divide(WrAer, rootdepth * 1000, out=np.zeros_like(WrAer), where=rootdepth!=0)

        # Calculate total available water and root zone depletion
        self.var.TAW = np.clip((WrFC - WrWP), 0, None)
        self.var.Dr = np.clip((WrFC - Wr), 0, None)
        self.var.Wr = np.copy(Wr)

