#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class TemperatureStress(object):
    def __init__(self, TemperatureStress_variable):
        self.var = TemperatureStress_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.Kst_Bio = np.copy(arr_zeros)
        self.var.Kst_PolH = np.copy(arr_zeros)
        self.var.Kst_PolC = np.copy(arr_zeros)

    def temperature_stress_biomass(self):
        """Function to calculate temperature stress coefficient 
        affecting biomass growth
        """
        KsBio_up = 1
        KsBio_lo = 0.02
        fshapeb = -1 * (np.log(((KsBio_lo * KsBio_up) - 0.98 * KsBio_lo) / (0.98 * (KsBio_up - KsBio_lo))))
        cond1 = (self.var.BioTempStress == 0)
        self.var.Kst_Bio[cond1] = 1
        cond2 = (self.var.BioTempStress == 1)
        cond21 = (cond2 & (self.var.GDD >= self.var.GDD_up))
        self.var.Kst_Bio[cond21] = 1
        cond22 = (cond2 & (self.var.GDD <= self.var.GDD_lo))
        self.var.Kst_Bio[cond22] = 0
        cond23 = (cond2 & np.logical_not(cond21 | cond22))
        GDDrel_divd = (self.var.GDD - self.var.GDD_lo)
        GDDrel_divs = (self.var.GDD_up - self.var.GDD_lo)
        GDDrel = np.divide(GDDrel_divd, GDDrel_divs, out=np.zeros_like(GDDrel_divs), where=GDDrel_divs!=0)
        Kst_Bio_divd = (KsBio_up * KsBio_lo)
        Kst_Bio_divs = (KsBio_lo + (KsBio_up - KsBio_lo) * np.exp(-fshapeb * GDDrel))
        self.var.Kst_Bio[cond23] = np.divide(Kst_Bio_divd, Kst_Bio_divs, out=np.zeros_like(Kst_Bio_divs), where=Kst_Bio_divs!=0)[cond23]
        self.var.Kst_Bio[cond23] = (self.var.Kst_Bio - KsBio_lo * (1 - GDDrel))[cond23]

    def temperature_stress_heat(self, KsPol_up, KsPol_lo):
        """Function to calculate effects of heat stress on 
        pollination
        """
        tmax = self.var.tmax[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
        cond3 = (self.var.PolHeatStress == 0)
        self.var.Kst_PolH[cond3] = 1
        cond4 = (self.var.PolHeatStress == 1)
        cond41 = (cond4 & (tmax <= self.var.Tmax_lo))
        self.var.Kst_PolH[cond41] = 1
        cond42 = (cond4 & (tmax >= self.var.Tmax_up))
        self.var.Kst_PolH[cond42] = 0
        cond43 = (cond4 & np.logical_not(cond41 | cond42))
        Trel_divd = (tmax - self.var.Tmax_lo)
        Trel_divs = (self.var.Tmax_up - self.var.Tmax_lo)
        Trel = np.divide(Trel_divd, Trel_divs, out=np.zeros_like(Trel_divs), where=Trel_divs!=0)
        Kst_PolH_divd = (KsPol_up * KsPol_lo)
        Kst_PolH_divs = (KsPol_lo + (KsPol_up - KsPol_lo) * np.exp(-self.var.fshape_b * (1 - Trel)))
        self.var.Kst_PolH[cond43] = np.divide(Kst_PolH_divd, Kst_PolH_divs, out=np.zeros_like(Kst_PolH_divs), where=Kst_PolH_divs!=0)[cond43]
        # return Kst_PolH
        

    def temperature_stress_cold(self, KsPol_up, KsPol_lo):
        """Function to calculate effects of cold stress on 
        pollination
        """
        tmin = self.var.tmin[None,:,:] * np.ones((self.var.nCrop))[:,None,None]
        cond5 = (self.var.PolColdStress == 0)
        self.var.Kst_PolC[cond5] = 1
        cond6 = (self.var.PolColdStress == 1)
        cond61 = (cond6 & (tmin >= self.var.Tmin_up))
        self.var.Kst_PolC[cond61] = 1
        cond62 = (cond6 & (tmin <= self.var.Tmin_lo))
        self.var.Kst_PolC[cond62] = 0
        Trel_divd = (self.var.Tmin_up - tmin)
        Trel_divs = (self.var.Tmin_up - self.var.Tmin_lo)
        Trel = np.divide(Trel_divd, Trel_divs, out=np.zeros_like(Trel_divs), where=Trel_divs!=0)
        Kst_PolC_divd = (KsPol_up * KsPol_lo)
        Kst_PolC_divs = (KsPol_lo + (KsPol_up - KsPol_lo) * np.exp(-self.var.fshape_b * (1 - Trel)))
        self.var.Kst_PolC[cond62] = np.divide(Kst_PolC_divd, Kst_PolC_divs, out=np.zeros_like(Kst_PolC_divs), where=Kst_PolC_divs!=0)[cond62]
        
    def dynamic(self):
        """Function to calculate temperature stress coefficients"""

        self.temperature_stress_biomass()

        KsPol_up = 1
        KsPol_lo = 0.001
        self.temperature_stress_heat(KsPol_up, KsPol_lo)
        self.temperature_stress_cold(KsPol_up, KsPol_lo)
