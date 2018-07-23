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
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLon, self.var.nLat))
        self.var.thRZ_Act = np.copy(arr_zeros)
        self.var.thRZ_Sat = np.copy(arr_zeros)
        self.var.thRZ_Fc = np.copy(arr_zeros)
        self.var.thRZ_Wp = np.copy(arr_zeros)
        self.var.thRZ_Dry = np.copy(arr_zeros)
        # self.thRZ_Aer = np.copy(arr_zeros)
        self.var.TAW = np.copy(arr_zeros)
        self.var.Dr = np.copy(arr_zeros)
        self.var.Wr = np.copy(arr_zeros)

    def fraction_of_compartment_in_root_zone(self):
        # Calculate root zone water content and available water
        rootdepth = np.maximum(self.var.Zmin, self.var.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100
        rootdepth_comp = np.broadcast_to(rootdepth[:,None,:,:], self.var.th.shape)
        comp_sto = (np.round((self.var.dzsum_xy - self.var.dz_xy) * 1000) < np.round(rootdepth_comp * 1000))

        # Fraction of compartment covered by root zone (zero in compartments
        # NOT covered by the root zone)
        factor = 1 - ((self.var.dzsum_xy - rootdepth_comp) / self.var.dz_xy)
        factor = np.clip(factor, 0, 1)
        factor[np.logical_not(comp_sto)] = 0

        self.var.RootDepth = rootdepth
        self.var.RootFact = factor
        # return rootdepth, factor

    def dynamic(self):
        """Function to calculate actual and total available water in the 
        root zone at current time step
        """
        self.fraction_of_compartment_in_root_zone()

        # Water storages in root zone (mm) - initially compute value in each
        # compartment, then sum to get overall root zone storages
        Wr_comp = self.var.RootFact * 1000 * self.var.th * self.var.dz_xy
        WrS_comp = self.var.RootFact * 1000 * self.var.th_s_comp * self.var.dz_xy
        WrFC_comp = self.var.RootFact * 1000 * self.var.th_fc_comp * self.var.dz_xy
        WrWP_comp = self.var.RootFact * 1000 * self.var.th_wp_comp * self.var.dz_xy
        WrDry_comp = self.var.RootFact * 1000 * self.var.th_dry_comp * self.var.dz_xy
        
        Wr = np.sum(Wr_comp, axis=1)
        Wr[Wr < 0] = 0
        WrS = np.sum(WrS_comp, axis=1)
        WrFC = np.sum(WrFC_comp, axis=1)
        WrWP = np.sum(WrWP_comp, axis=1)
        WrDry = np.sum(WrDry_comp, axis=1)

        # Convert depths to m3/m3
        self.var.thRZ_Act = np.divide(Wr, self.var.RootDepth * 1000, out=np.zeros_like(Wr), where=self.var.RootDepth!=0)
        self.var.thRZ_Sat = np.divide(WrS, self.var.RootDepth * 1000, out=np.zeros_like(WrS), where=self.var.RootDepth!=0)
        self.var.thRZ_Fc  = np.divide(WrFC, self.var.RootDepth * 1000, out=np.zeros_like(WrFC), where=self.var.RootDepth!=0)
        self.var.thRZ_Wp  = np.divide(WrWP, self.var.RootDepth * 1000, out=np.zeros_like(WrWP), where=self.var.RootDepth!=0)
        self.var.thRZ_Dry = np.divide(WrDry, self.var.RootDepth * 1000, out=np.zeros_like(WrDry), where=self.var.RootDepth!=0)

        # # Water storage in root zone at aeration stress threshold (mm)
        # WrAer_comp = self._factor * 1000 * (self.var.th_s_comp - (self.var.Aer / 100)) * dz
        # WrAer = np.sum(WrAer_comp, axis=0)
        # self.var.thRZ_Aer = np.divide(WrAer, self._rootdepth * 1000, out=np.zeros_like(WrAer), where=self._rootdepth!=0)

        # Calculate total available water and root zone depletion
        self.var.TAW = np.clip((WrFC - WrWP), 0, None)
        self.var.Dr = np.clip((WrFC - Wr), 0, None)
        self.var.Wr = np.copy(Wr)

class AQRootZoneWater(RootZoneWater):
    def initial(self):
        super(AQRootZoneWater, self).initial()
        self.var.thRZ_Aer = np.zeros_like(self.var.thRZ_Act)

    def dynamic(self):
        super(AQRootZoneWater, self).dynamic()
        WrAer_comp = self.var.RootFact * 1000 * (self.var.th_s_comp - (self.var.Aer / 100)) * self.var.dz_xy
        WrAer = np.sum(WrAer_comp, axis=1)
        self.var.thRZ_Aer = np.divide(WrAer, self.var.RootDepth * 1000, out=np.zeros_like(WrAer), where=self.var.RootDepth!=0)
        
        
class FAO56RootZoneWater(RootZoneWater):
    """Class to represent root zone water in FAO56 model"""

