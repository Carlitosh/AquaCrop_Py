#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import numpy as np

import logging
logger = logging.getLogger(__name__)

class Irrigation(object):
    """Class to represent irrigation activities"""

    def __init__(self, Irrigation_variable):
        self.var = Irrigation_variable

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.Irr = np.copy(arr_zeros)
        self.var.IrrCum = np.copy(arr_zeros)
        self.var.IrrNetCum = np.copy(arr_zeros)

    def reset_initial_conditions(self):
        self.var.IrrCum[self.var.GrowingSeasonDayOne] = 0
        self.var.IrrNetCum[self.var.GrowingSeasonDayOne] = 0
        
    def dynamic(self):
        """Function to get irrigation depth for the current day"""
        # reset initial conditions
        if np.any(self.var.GrowingSeasonDayOne):
            self.reset_initial_conditions()
        
        SMT = np.concatenate((self.var.SMT1[None,:],
                              self.var.SMT2[None,:],
                              self.var.SMT3[None,:],
                              self.var.SMT4[None,:]), axis=0)

        # Determine adjustment for inflows and outflows on current day
        cond1 = (self.var.thRZ_Act > self.var.thRZ_Fc)
        rootdepth = np.maximum(self.var.Zmin, self.var.Zroot)
        AbvFc = ((self.var.thRZ_Act - self.var.thRZ_Fc) * 1000 * rootdepth)
        AbvFc[np.logical_not(cond1)] = 0
        WCadj = self.var.ETpot - self.var.precipitation + self.var.Runoff - AbvFc
        # WCadj = self.var.Tpot + self.var.Epot - self.var.precipitation + self.var.Runoff - AbvFc

        # Determine irrigation depth (mm/day) to be applied
        cond2 = (self.var.GrowingSeasonIndex & (self.var.IrrMethod == 0))
        self.var.Irr[cond2] = 0
        
        # If irrigation is based on soil moisture, get the soil moisture
        # target for the current growth stage and determine threshold to
        # initiate irrigation
        cond3 = (self.var.GrowingSeasonIndex & np.logical_not(cond2) & (self.var.IrrMethod == 1))
        I,J,K = np.ogrid[:self.var.nCrop,:self.var.nLat,:self.var.nLon]
        growth_stage_index = self.var.GrowthStage.astype(int) - 1
        SMT = SMT[growth_stage_index,I,J,K]
        IrrThr = np.round(((1 - SMT / 100) * self.var.TAW), 3)

        # If irrigation is based on a fixed interval, get number of days in
        # growing season so far (subtract 1 so that we always irrigate first
        # on day 1 of each growing season)
        cond4 = (self.var.GrowingSeasonIndex & (self.var.IrrMethod == 2))
        nDays = self.var.DAP - 1

        # Working on a copy, adjust depletion for inflows and outflows - same
        # for both soil moisture and interval based irrigation
        cond5 = (cond3 | cond4)
        Dr = np.copy(self.var.Dr)
        Dr[cond5] += WCadj[cond5]
        Dr[cond5] = np.clip(Dr, 0, None)[cond5]
        cond6 = ((cond3 & (Dr > IrrThr)) | (cond4 & ((nDays % self.var.IrrInterval) == 0)))        
        IrrReq = np.copy(Dr)

        EffAdj = ((100 - self.var.AppEff) + 100) / 100
        IrrReq *= EffAdj
        self.var.Irr[cond6] = np.clip(IrrReq, 0, self.var.MaxIrr)[cond6]
        cond7 = (cond5 & np.logical_not(cond6))
        self.var.Irr[cond7] = 0

        # If irrigation is based on a pre-defined schedule then the irrigation
        # requirement for each crop is read from a netCDF file. Note that if
        # the option 'irrScheduleFileNC' is None, then nothing will be imported
        # and the irrigation requirement will be zero
        cond8 = (self.var.GrowingSeasonIndex & (self.var.IrrMethod == 3))
        if self.var.irrScheduleFileNC != None:
            IrrReq = vos.netcdf2PCRobjClone(self.var.irrScheduleFileNC,\
                                            "irrigation_schedule",\
                                            str(currTimeStep.fulldate),\
                                            useDoy = method_for_time_index,\
                                            cloneMapFileName = self.var.cloneMap,\
                                            LatitudeLongitude = True)
            self.var.Irr[cond8] = IrrReq[cond8]

        # Note that if irrigation is based on net irrigation then it is
        # performed after calculation of transpiration and hence is set to zero
        # at this point (not strictly necessary because Irr is initialized to
        # zero, but included for completeness)
        cond9 = (self.var.GrowingSeasonIndex & (self.var.IrrMethod == 4))
        self.var.Irr[cond9] = 0
        self.var.Irr[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        self.var.IrrCum += self.var.Irr
