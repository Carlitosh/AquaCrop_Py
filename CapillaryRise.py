#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

# NB: this class has not been tested!

import numpy as np

import logging
logger = logging.getLogger(__name__)

class CapillaryRise(object):
    def __init__(self, CapillaryRise_variable):
        self.var = CapillaryRise_variable

    def initial(self):
        self.var.CrTot = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

    def maximum_capillary_rise(ksat, aCR, bCR, zGW, z):
        """Function to compute maximum capillary rise for a given soil 
        profile

        Args:
          ksat : saturated hydraulic conductivity of the soil layer
          aCR  : ... of the soil layer
          bCR  : ... of the soil layer
          zGW  : depth of groundwater table
          z    : depth to midpoint of the soil layer

        """
        MaxCR = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        cond1 = ((ksat > 0) & (zGW > 0) & ((zGW - z) < 4))
        cond11 = (cond1 & (z >= zGW))
        MaxCR[cond11] = 99            
        cond12 = (cond1 & np.logical_not(cond11))
        MaxCR_log = np.log(zGW - z, out=np.zeros((MaxCR.shape)), where=cond12)
        MaxCR[cond12] = (np.exp((MaxCR_log - bCR) / aCR))[cond12]
        MaxCR = np.clip(MaxCR, None, 99)
        return MaxCR

    def store_water_from_capillary_rise(th, th_fc, th_fc_adj, th_wp, fshape_cr, MaxCR, flux_out, zGW, zBot, dz):

        thnew = np.copy(th)

        cond1 = ((np.round(MaxCR * 1000) > 0) & (np.round(flux_out * 1000) == 0))

        # calculate driving force
        Df = driving_force(th, th_fc_adj, th_wp, fshape_cr)

        # Calculate relative hydraulic conductivity
        Krel = relative_hydraulic_conductivity(th, th_fc, th_wp)

        # Check if room is available to store water from capillary rise
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        dth = np.copy(arr_zeros)
        dth[cond1] = (th_fc_adj - th)[cond1]                
        dth = np.round((dth * 1000) / 1000)

        # Store water if room is available
        dthMax = np.copy(arr_zeros)
        CRcomp = np.copy(arr_zeros)
        cond15 = (cond1 & (dth > 0) & ((zBot - dz / 2) < zGW))

        dthMax[cond35] = (Krel * Df * MaxCR / (1000 * dz))[cond35]
        cond151 = (cond15 & (dth >= dthMax))
        thnew[cond151] += dthMax[cond151]
        CRcomp[cond151] = (dthMax * 1000 * dz)[cond151]
        MaxCR[cond151] = 0

        cond152 = (cond35 & np.logical_not(cond151))
        thnew[cond152] = th_fc_adj[cond152]
        CRcomp[cond152] = (dth * 1000 * dz)[cond152]
        MaxCR[cond152] = ((Krel * MaxCR) - CRcomp)[cond152]    
        return thnew, CRcomp, MaxCR
    
    def dynamic(self):
        """Function to calculate capillary rise from a shallow 
        groundwater table
        """
        if self.var.WaterTable:
            
            # dims = self.var.th.shape
            # nc, nr, nlat, nlon = dims[0], dims[1], dims[2], dims[3]
            zBot = np.sum(self.var.dz)
            zBotMid = zBot - (self.var.dz[-1] / 2)  # depth to midpoint of bottom layer
            thnew = np.copy(self.var.th)

            # Get maximum capillary rise for bottom compartment
            MaxCR = self.maximum_capillary_rise(
                self.var.ksat[:,-1,...],
                self.var.aCR[:,-1,...],
                self.var.bCR[:,-1,...],
                self.var.zGW,
                zBotMid)

            # Check for restrictions on upward flow caused by properties of
            # compartments that are not modelled in the soil water balance

            # Find top of next soil layer that is not within modelled soil profile
            zTopLayer = np.cumsum(self.var.zLayer[np.arange(0, (self.var.layerIndex[-1] + 1))])
            layeri = self.var.layerIndex[-1]  # layer number of bottom compartment
            LimCR = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

            # TODO: should this be nlayer, rather than nlayer-1
            while np.any(zTopLayer < self.var.zGW) & (layeri < (self.var.nLayer - 1)):
                layeri += 1
                LimCR = self.maximum_capillary_rise(
                    self.var.ksat[:,layeri,...],
                    self.var.aCR[:,layeri,...],
                    self.var.bCR[:,layeri,...],
                    self.var.zGW,
                    zTopLayer)
                MaxCR = np.clip(MaxCR, None, LimCR)
                zTopLayer += self.var.zLayer[layeri]

            compi = self.var.nComp - 1
            CrTot = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
            while ((np.any(np.round(MaxCR * 1000) > 0)) & (np.any(np.round(self.var.FluxOut[:,compi,...] * 1000) == 0)) & (compi >= 0)):

                cond1 = ((np.round(MaxCR * 1000) > 0) & (np.round(self.var.FluxOut[:,compi,...] * 1000) == 0))
                thnew_comp, CRcomp, MaxCR = self.store_water_from_capillary_rise(
                    self.var.th[:,compi,...],
                    self.var.th_fc_comp[:,compi,...],
                    self.var.th_fc_adj[:,compi,...],
                    self.var.th_wp_comp[:,compi,...],
                    self.var.fshape_cr_comp[:,compi,...],
                    MaxCR,
                    self.var.FluxOut[:,compi,...],
                    self.var.zGW,
                    zBot,
                    self.var.dz[compi])

                thnew[:,compi,...][cond1] = thnew_comp[cond1]
                CrTot[cond1] += CRcomp[cond1]

                # Update bottom elevation of compartment and compartment counter
                zBot -= self.var.dz[compi]
                compi -= 1

                # Update restriction on maximum capillary rise
                if compi >= 0:
                    zBotMid = zBot - (self.var.dz[compi] / 2)
                    LimCR = self.maximum_capillary_rise(
                        self.var.ksat_comp[:,compi,...],
                        self.var.aCR_comp[:,compi,...],
                        self.var.bCR_comp[:,compi,...],
                        self.var.zGW,
                        zBotMid)
                    cond11 = (cond1 & (MaxCR > LimCR))
                    MaxCR[cond11] = LimCR[cond11]

            self.var.th = np.copy(thnew)
            self.var.CrTot = CrTot
