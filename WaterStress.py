#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class WaterStress(object):
    def __init__(self, WaterStress_variable):
        self.var = WaterStress_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
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

        et0 = (self.var.referencePotET[None,:,:] * np.ones((self.var.nCrop))[:,None,None])
        
        dims = et0.shape
        nr, nlat, nlon = dims[0], dims[1], dims[2]  # TODO: get rid of this

        # Adjust stress thresholds for Et0 on current day (don't do this for
        # pollination water stress coefficient)
        cond1 = (self.var.ETadj == 1)
        for stress in range(3):
            p_up[stress,:][cond1] = (p_up[stress,:] + (0.04 * (5 - et0)) * (np.log10(10 - 9 * p_up[stress,:])))[cond1]
            p_lo[stress,:][cond1] = (p_lo[stress,:] + (0.04 * (5 - et0)) * (np.log10(10 - 9 * p_lo[stress,:])))[cond1]

        # Adjust senescence threshold if early senescence triggered
        if beta:
            cond2 = (self.var.tEarlySen > 0)
            p_up[2,:][cond2] = (p_up[2,:] * (1. - (self.var.beta / 100.)))[cond2]

        # Limit adjusted values
        p_up = np.clip(p_up, 0, 1)
        p_lo = np.clip(p_lo, 0, 1)

        # Calculate relative depletion
        Drel = np.zeros((4, nr, nlat, nlon))
        # No water stress
        cond1 = (self.var.Dr <= (p_up * self.var.TAW))
        Drel[cond1] = 0

        # Partial water stress
        cond2 = (self.var.Dr > (p_up * self.var.TAW)) & (self.var.Dr < (p_lo * self.var.TAW)) & np.logical_not(cond1)
        x1 = p_lo - np.divide(self.var.Dr, self.var.TAW, out=np.zeros_like(Drel), where=self.var.TAW!=0)
        x2 = p_lo - p_up
        Drel[cond2] = (1 - np.divide(x1, x2, out=np.zeros_like(Drel), where=x2!=0))[cond2]

        # Full water stress
        cond3 = (self.var.Dr >= (p_lo * self.var.TAW)) & np.logical_not(cond1 | cond2)
        Drel[cond3] = 1         

        # Calculate root zone stress coefficients
        idx = np.arange(0,3)
        x1 = np.exp(Drel[idx,:] * fshape_w[idx,:]) - 1.
        x2 = np.exp(fshape_w[idx,:]) - 1.
        Ks = (1. - np.divide(x1, x2, out=np.zeros_like(x2), where=x2!=0))
        # print 'Ksw_Sen in WaterStress %f' % Ks[2,0,0,0]
        # Water stress coefficients (leaf expansion, stomatal closure,
        # senescence, pollination failure)
        self.var.Ksw_Exp = np.copy(Ks[0,:])
        self.var.Ksw_Sto = np.copy(Ks[1,:])
        self.var.Ksw_Sen = np.copy(Ks[2,:])
        self.var.Ksw_Pol = 1 - Drel[3,:]

        # Mean water stress coefficient for stomatal closure
        self.var.Ksw_StoLin = 1 - Drel[1,:]
