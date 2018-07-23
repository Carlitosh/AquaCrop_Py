#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class RainfallPartition(object):
    """Class to infiltrate incoming water"""
    
    def __init__(self, RainfallPartition_variable):
        self.var = RainfallPartition_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))        
        self.var.Runoff = np.copy(arr_zeros)
        self.var.Infl = np.copy(arr_zeros)
        
    def dynamic(self):
        """Function to partition rainfall into surface runoff and 
        infiltration using the curve number approach.
        """
        # Add crop dimension to precipitation
        P = self.var.precipitation[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
        zcn = self.var.zCN[:,None,:,:] * np.ones((self.var.nComp))[None,:,None,None]

        cond1 = ((self.var.Bunds == 0) | (self.var.zBund < 0.001))

        # Adjust curve number
        cond11 = (cond1 & (self.var.AdjCN == 1))

        # Check which compartment cover depth of top soil used to adjust
        # curve number
        comp_sto = (np.round(self.var.dzsum_xy * 1000) <= np.round(zcn * 1000))
        cond111 = np.all((comp_sto == False), axis=1)
        cond111 = np.broadcast_to(cond111[:,None,:,:], comp_sto.shape)
        comp_sto[cond111] = True

        # Calulcate weighting factors by compartment
        dzsum = np.copy(self.var.dzsum_xy)
        dzsum[dzsum > zcn] = zcn[dzsum > zcn]
        wx = (1.016 * (1 - np.exp(-4.16 * (dzsum / zcn))))

        # xx is wx for the overlying layer, with the top layer equal to zero
        xx = np.concatenate((np.zeros((self.var.nCrop, 1, self.var.nLat, self.var.nLon)), wx[:,:-1,...]), axis=1)
        wrel = np.clip((wx - xx), 0, 1)

        # Calculate relative wetness of topsoil
        thnew = np.maximum(self.var.th_wp_comp, self.var.th)

        # Multiply by comp_sto to ensure that compartments not used for
        # curve number adjustment are set to zero
        wet_top_comp = wrel * ((self.var.th - self.var.th_wp_comp) / (self.var.th_fc_comp - self.var.th_wp_comp)) * comp_sto
        wet_top = np.sum(wet_top_comp, axis=1)  # sum along compartment dimension
        wet_top = np.clip(wet_top, 0, 1)

        # Make adjustment to curve number
        CN = np.copy(self.var.CN)
        CN[cond11] = np.round(self.var.CNbot + (self.var.CNtop - self.var.CNbot) * wet_top)[cond11]

        # Partition rainfall into runoff and infiltration
        S = (25400. / CN) - 254
        term = (P - ((5. / 100) * S))

        cond12 = (cond1 & (term <= 0))
        self.var.Runoff[cond12] = 0
        self.var.Infl[cond12] = P[cond12]

        cond13 = (cond1 & np.logical_not(cond12))
        self.var.Runoff[cond13] = ((term ** 2) / (P + (1 - (5. / 100)) * S))[cond13]
        self.var.Infl[cond13] = (P - self.var.Runoff)[cond13]

        # If bunds are present then there is no runoff
        cond2 = np.logical_not(cond1)
        self.var.Runoff[cond2] = 0
        self.var.Infl[cond2] = P[cond2]
