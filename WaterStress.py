#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

from crop_growth_funs import *

import logging
logger = logging.getLogger(__name__)

class WaterStress(object):
    def __init__(self, WaterStress_variable):
        self.var = WaterStress_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.var.Ksw_Exp = np.copy(arr_zeros)
        self.var.Ksw_Sto = np.copy(arr_zeros)
        self.var.Ksw_Sen = np.copy(arr_zeros)
        self.var.Ksw_Pol = np.copy(arr_zeros)
        self.var.Ksw_StoLin = np.copy(arr_zeros)
        
    def dynamic(self, beta):
        """Function to calculate water stress coefficients"""        
        p_up = np.concatenate((self.var.p_up1[None,:], self.var.p_up2[None,:], self.var.p_up3[None,:], self.var.p_up4[None,:]), axis=0)
        p_lo = np.concatenate((self.var.p_lo1[None,:], self.var.p_lo2[None,:], self.var.p_lo3[None,:], self.var.p_lo4[None,:]), axis=0)
        fshape_w = np.concatenate((self.var.fshape_w1[None,:], self.var.fshape_w2[None,:], self.var.fshape_w3[None,:], self.var.fshape_w4[None,:]), axis=0)

        et0 = (self.var.referencePotET[None,:,:] * np.ones((self.var.nRotation))[:,None,None])
        self.var.Ksw_Exp, self.var.Ksw_Sto, self.var.Ksw_Sen, self.var.Ksw_Pol, self.var.Ksw_StoLin = water_stress(
            p_lo, p_up, fshape_w, et0, self.var.ETadj, self.var.tEarlySen,
            beta, self.var.beta, self.var.TAW, self.var.Dr)
