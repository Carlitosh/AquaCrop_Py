#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class GrowingSeason(object):
    """Class to represent crop growing season"""
    def __init__(self, GrowingSeason_variable):
        self.var = GrowingSeason_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        self.var.GrowingSeasonDayOne = np.copy(arr_zeros.astype(bool))
        self.var.GrowingSeasonIndex = np.copy(arr_zeros.astype(bool))
        self.var.DAP = np.copy(arr_zeros)

    def reset_initial_conditions(self):
        self.var.DAP[self.var.GrowingSeasonDayOne] = 0

    def dynamic(self):
        self.var.GrowingSeasonDayOne = self.var._modelTime.doy == self.var.PlantingDate
        self.reset_initial_conditions()
        self.var.GrowingSeasonIndex *= np.logical_not(self.var.CropDead | self.var.CropMature)
        self.var.GrowingSeasonIndex[self.var.GrowingSeasonDayOne] = True
        self.var.DAP[self.var.GrowingSeasonIndex] += 1
        self.var.DAP[np.logical_not(self.var.GrowingSeasonIndex)] = 0
