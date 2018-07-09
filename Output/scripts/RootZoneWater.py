#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from hydrology_funs import *

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
        thRZ, TAW, Dr, Wr = root_zone_water(
            self.var.th,
            self.var.th_s_comp,
            self.var.th_fc_comp,
            self.var.th_wp_comp,
            self.var.th_dry_comp,
            self.var.dz,
            self.var.Zmin,
            self.var.Zroot,
            self.var.Aer)

        self.var.thRZ_Act = thRZ['Act']
        self.var.thRZ_Sat = thRZ['Sat']
        self.var.thRZ_Fc  = thRZ['Fc']
        self.var.thRZ_Wp  = thRZ['Wp']
        self.var.thRZ_Dry = thRZ['Dry']
        self.var.thRZ_Aer = thRZ['Aer']
        self.var.TAW = TAW
        self.var.Dr = Dr
        self.var.Wr = Wr
