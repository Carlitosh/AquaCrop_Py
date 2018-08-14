#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# AquaCrop
#
import os
import numpy as np
import pcraster as pcr
import VirtualOS as vos
import netCDF4 as nc
import datetime as datetime
import calendar as calendar

class CropParameters(object):
    
    def __init__(self, CropParameters_variable):
        self.var = CropParameters_variable
        self.var.nCrop = int(self.var._configuration.cropOptions['nCrop'])
        self.var.cropParameterFileNC = str(self.var._configuration.cropOptions['cropParameterNC'])
        self.var.crop_parameters_to_read = []
        self.var.crop_parameters_to_compute = []

    def initial(self):
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.GrowingSeasonIndex = np.copy(arr_zeros.astype(bool))
        self.var.GrowingSeasonDayOne = np.copy(arr_zeros.astype(bool))
        self.var.DAP = np.copy(arr_zeros)
        
    def read(self):
        """Function to read crop input parameters"""
        if len(self.var.crop_parameters_to_read) > 0:
            for param in self.var.crop_parameters_to_read:
                # nm = '_' + param
                vars(self.var)[param] = vos.netcdf2PCRobjCloneWithoutTime(
                    self.var.cropParameterFileNC,
                    param,
                    cloneMapFileName=self.var.cloneMap)
        
    def adjust_planting_and_harvesting_date(self):

        if self.var._modelTime.timeStepPCR == 1 or self.var._modelTime.doy == 1:
            
            pd = np.copy(self.var.PlantingDate)
            hd = np.copy(self.var.HarvestDate)
            st = self.var._modelTime.currTime
            sd = self.var._modelTime.currTime.timetuple().tm_yday

            # adjust values for leap year (objective is to preserve date)
            isLeapYear1 = calendar.isleap(st.year)
            # isLeapYear2 = calendar.isleap(st.year + 1)
            pd[(isLeapYear1 & (pd >= 60))] += 1  # TODO: check these
            hd[(isLeapYear1 & (hd >= 60) & (hd > pd))] += 1
            # hd[(isLeapYear2 & (hd >= 60) & (hd < pd))] += 1
            
            self.var.PlantingDateAdj = np.copy(pd)
            self.var.HarvestDateAdj = np.copy(hd)

    def update_growing_season(self):

        # TODO: check both PlantingDateAdj and HarvestDateAdj are not the same!        
        cond1 = ((self.var.PlantingDateAdj < self.var.HarvestDateAdj)
                 & ((self.var.PlantingDateAdj <= self.var._modelTime.doy)
                    & (self.var._modelTime.doy <= self.var.HarvestDateAdj)))
        cond2 = ((self.var.PlantingDateAdj > self.var.HarvestDateAdj)
                 & ((self.var.PlantingDateAdj <= self.var._modelTime.doy)
                    | (self.var._modelTime.doy <= self.var.HarvestDateAdj)))

        # ***TODO: introduce another condition to check whether current growing season starts before simulation start time***
        
        hd_arr = (np.datetime64(str(datetime.datetime(self.var._modelTime.year, 1, 1))) + np.array(self.var.HarvestDateAdj - 1, dtype='timedelta64[D]'))        
        cond3 = hd_arr <= np.datetime64(str(self.var._modelTime.endTime))
        
        self.var.GrowingSeason = ((cond1 & cond3) | (cond2 & cond3))

        # print self.var.PlantingDateAdj[:,10,10]
        # print self.var.HarvestDateAdj[:,10,10]
                
        self.var.GrowingSeasonIndex = np.copy(self.var.GrowingSeason)
        self.var.GrowingSeasonIndex *= np.logical_not(self.var.CropDead | self.var.CropMature)

        self.var.GrowingSeasonDayOne = self.var._modelTime.doy == self.var.PlantingDateAdj
        # self.var.GrowingSeasonDayOne = self.var._modelTime.doy == self.var.PlantingDate
        self.var.DAP[self.var.GrowingSeasonDayOne] = 0
        self.var.GrowingSeasonIndex[self.var.GrowingSeasonDayOne] = True

        self.var.DAP[self.var.GrowingSeasonIndex] += 1
        self.var.DAP[np.logical_not(self.var.GrowingSeasonIndex)] = 0
        
class AQCropParameters(CropParameters):

    def initial(self):

        super(AQCropParameters, self).initial()

        # TODO: move these to main model class?
        self.var.CalendarType = int(self.var._configuration.cropOptions['CalendarType'])
        self.var.SwitchGDD = bool(int(self.var._configuration.cropOptions['SwitchGDD']))
        self.var.GDDmethod = int(self.var._configuration.cropOptions['GDDmethod'])
        
        # Declare variables
        self.var.crop_parameters_to_read = [
            'CropType','PlantingDate','HarvestDate','Emergence','MaxRooting',
            'Senescence','Maturity','HIstart','Flowering','YldForm',
            'PolHeatStress','PolColdStress','BioTempStress','PlantPop',
            'Determinant','ETadj','LagAer','Tbase','Tupp','Tmax_up','Tmax_lo',
            'Tmin_up','Tmin_lo','GDD_up','GDD_lo','fshape_b','PctZmin','Zmin',
            'Zmax','fshape_r','fshape_ex','SxTopQ','SxBotQ','a_Tr','SeedSize',
            'CCmin','CCx','CDC','CGC','Kcb','fage','WP','WPy','fsink','bsted',
            'bface','HI0','HIini','dHI_pre','a_HI','b_HI','dHI0','exc',
            'MaxFlowPct','p_up1','p_up2','p_up3','p_up4','p_lo1','p_lo2','p_lo3',
            'p_lo4','fshape_w1','fshape_w2','fshape_w3','fshape_w4','Aer','beta',
            'GermThr']

        self.var.crop_parameters_to_compute = [
            'CC0','SxTop','SxBot','fCO2','tLinSwitch','dHILinear','HIGC',
            'CanopyDevEnd','CanopyDevEndCD','Canopy10Pct','Canopy10PctCD','MaxCanopy',
            'MaxCanopyCD','HIstartCD','HIend','HIendCD','YldFormCD',
            'FloweringEnd','Flowering','FloweringCD','CGC','CDC',
            'PlantingDateAdj','HarvestDateAdj',
            'CurrentConc']  # TODO: CGC and CDC are in both list - try removing them here?

        self.var.crop_parameter_names = self.var.crop_parameters_to_read + self.var.crop_parameters_to_compute
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        for param in self.var.crop_parameter_names:
            vars(self.var)[param] = np.copy(arr_zeros)
            
        self.read()
        self.compute_crop_parameters()

    def compute_crop_parameters(self):
        """Function to compute additional crop variables required to run 
        AquaCrop"""
        
        # The following adapted from lines 160-240 of AOS_ComputeVariables.m

        # Fractional canopy cover size at emergence
        self.var.CC0 = np.round(10000. * (self.var.PlantPop * self.var.SeedSize) * 10 ** -8) / 10000
        
        # Root extraction terms
        S1 = np.copy(self.var.SxTopQ)
        S2 = np.copy(self.var.SxBotQ)
        self.var.SxTop = np.zeros_like(self.var.SxTopQ)
        self.var.SxBot = np.zeros_like(self.var.SxBotQ)

        cond1 = (S1 == S2)
        self.var.SxTop[cond1] = S1[cond1]
        self.var.SxBot[cond1] = S2[cond1]

        cond2 = np.logical_not(cond1)
        cond21 = (cond2 & (self.var.SxTopQ < self.var.SxBotQ))
        S1[cond21] = self.var.SxBotQ[cond21]
        S2[cond21] = self.var.SxTopQ[cond21]
        xx = 3 * (S2 / (S1 - S2))
        SS1 = np.zeros_like(S1)
        SS2 = np.zeros_like(S2)
        # SS1 = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        # SS2 = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

        cond22 = (cond2 & (xx < 0.5))
        SS1[cond22] = ((4 / 3.5) * S1)[cond22]
        SS2[cond22] = 0

        cond23 = (cond2 & np.logical_not(cond22))
        SS1[cond23] = ((xx + 3.5) * (S1 / (xx + 3)))[cond23]
        SS2[cond23] = ((xx - 0.5) * (S2 / xx))[cond23]

        cond24 = (cond2 & (self.var.SxTopQ > self.var.SxBotQ))
        self.var.SxTop[cond24] = SS1[cond24]
        self.var.SxBot[cond24] = SS2[cond24]

        cond25 = (cond2 & np.logical_not(cond24))
        self.var.SxTop[cond25] = SS2[cond25]
        self.var.SxBot[cond25] = SS1[cond25]

        # Crop calender
        self.compute_crop_calendar()

        # Harvest index growth coefficient
        self.calculate_HIGC()

        # Days to linear HI switch point
        self.calculate_HI_linear()

    def compute_water_productivity_adjustment_factor(self):
        """Function to calculate water productivity adjustment factor 
        for elevation in CO2 concentration"""

        # Get CO2 weighting factor
        fw = np.zeros_like(self.var.conc)
        cond1 = (self.var.conc > self.var.RefConc)
        cond11 = (cond1 & (self.var.conc >= 550))
        fw[cond11] = 1
        cond12 = (cond1 & np.logical_not(cond11))
        fw[cond12] = (1 - ((550 - self.var.conc) / (550 - self.var.RefConc)))[cond12]

        # Determine adjustment for each crop in first year of simulation
        fCO2 = ((self.var.conc / self.var.RefConc) /
                (1 + (self.var.conc - self.var.RefConc) * ((1 - fw)
                                           * self.var.bsted + fw
                                           * ((self.var.bsted * self.var.fsink)
                                              + (self.var.bface
                                                 * (1 - self.var.fsink))))))

        # Consider crop type
        ftype = (40 - self.var.WP) / (40 - 20)
        ftype = np.clip(ftype, 0, 1)
        fCO2 = 1 + ftype * (fCO2 - 1)
        
        self.var.fCO2[self.var.GrowingSeasonDayOne] = fCO2[self.var.GrowingSeasonDayOne]
        conc = (self.var.conc[None,:,:] * np.ones((self.var.nCrop))[:,None,None])
        self.var.CurrentConc[self.var.GrowingSeasonDayOne] = conc[self.var.GrowingSeasonDayOne]
        
    def calculate_HI_linear(self):
        """Function to calculate time to switch to linear harvest index 
        build-up, and associated linear rate of build-up. Only for 
        fruit/grain crops
        """
        # Determine linear switch point
        cond1 = (self.var.CropType == 3)
        ti = np.zeros_like(self.var.CropType)
        # ti = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        tmax = self.var.YldFormCD
        HIest = np.zeros_like(self.var.HIini)
        # HIest = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        HIprev = self.var.HIini

        # Iterate to find linear switch point
        while np.any((self.var.CropType == 3) & (HIest <= self.var.HI0) & (ti < tmax)):    
            ti += 1
            HInew = ((self.var.HIini * self.var.HI0) / (self.var.HIini + (self.var.HI0 - self.var.HIini) * np.exp(-self.var.HIGC * ti)))
            HIest = (HInew + (tmax - ti) * (HInew - HIprev))
            HIprev = HInew

        self.var.tLinSwitch = ti - 1  # Line 19 of AOS_CalculateHILinear.m
        # self.var.tLinSwitch[self.var.CropType != 3] = np.nan
            
        # Determine linear build-up rate
        cond1 = (self.var.tLinSwitch > 0)
        HIest[cond1] = ((self.var.HIini * self.var.HI0) / (self.var.HIini + (self.var.HI0 - self.var.HIini) * np.exp(-self.var.HIGC * self.var.tLinSwitch)))[cond1]
        HIest[np.logical_not(cond1)] = 0
        self.var.dHILinear = ((self.var.HI0 - HIest) / (tmax - self.var.tLinSwitch))  # dHILin will be set to nan in the same cells as tSwitch
                    
    def calculate_HIGC(self):
        """Function to calculate harvest index growth coefficient"""
        # Total yield formation days
        tHI = np.copy(self.var.YldFormCD)

        # Iteratively estimate HIGC
        self.var.HIGC = np.full((self.var.HIini.shape), 0.001)
        # self.var.HIGC = np.full((self.var.nCrop, self.var.nLat, self.var.nLon), 0.001)
        HIest = np.zeros_like(self.var.HIini)
        # HIest = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        while np.any(HIest < (0.98 * self.var.HI0)):
            cond1 = (HIest < (0.98 * self.var.HI0))
            self.var.HIGC += 0.001
            HIest = ((self.var.HIini * self.var.HI0) / (self.var.HIini + (self.var.HI0 - self.var.HIini) * np.exp(-self.var.HIGC * tHI)))

        self.var.HIGC[HIest >= self.var.HI0] -= 0.001

    def compute_crop_calendar(self):
       
        # "Time from sowing to end of vegetative growth period"
        cond1 = (self.var.Determinant == 1)
        self.var.CanopyDevEnd = np.copy(self.var.Senescence)
        self.var.CanopyDevEnd[cond1] = (np.round(self.var.HIstart + (self.var.Flowering / 2)))[cond1]

        # "Time from sowing to 10% canopy cover (non-stressed conditions)
        self.var.Canopy10Pct = np.round(self.var.Emergence + (np.log(0.1 / self.var.CC0) / self.var.CGC))

        # "Time from sowing to maximum canopy cover (non-stressed conditions)
        self.var.MaxCanopy = np.round(self.var.Emergence + (np.log((0.25 * self.var.CCx * self.var.CCx / self.var.CC0) / (self.var.CCx - (0.98 * self.var.CCx))) / self.var.CGC))

        # "Time from sowing to end of yield formation"
        self.var.HIend = self.var.HIstart + self.var.YldForm
        cond2 = (self.var.CropType == 3)

        # TODO: declare these in __init__
        arr_zeros = np.zeros_like(self.var.CropType)
        self.var.FloweringEnd = np.copy(arr_zeros)
        self.var.FloweringEndCD = np.copy(arr_zeros)
        self.var.FloweringCD = np.copy(arr_zeros)
        self.var.FloweringEnd[cond2] = (self.var.HIstart + self.var.Flowering)[cond2]
        self.var.FloweringEndCD[cond2] = self.var.FloweringEnd[cond2]
        self.var.FloweringCD[cond2] = self.var.Flowering[cond2]
        
        # Mode = self.var.CalendarType
        # if Mode == 1:
        if self.var.CalendarType == 1:
            # "Duplicate calendar values (needed to minimise if-statements when
            # switching between GDD and CD runs)
            self.var.EmergenceCD = np.copy(self.var.Emergence)      # only used in this function
            self.var.Canopy10PctCD = np.copy(self.var.Canopy10Pct)  # only used in this function
            self.var.MaxRootingCD = np.copy(self.var.MaxRooting)    # only used in this function
            self.var.SenescenceCD = np.copy(self.var.Senescence)    # only used in this function
            self.var.MaturityCD = np.copy(self.var.Maturity)        # only used in this function
            self.var.MaxCanopyCD = np.copy(self.var.MaxCanopy)
            self.var.CanopyDevEndCD = np.copy(self.var.CanopyDevEnd)
            self.var.HIstartCD = np.copy(self.var.HIstart)
            self.var.HIendCD = np.copy(self.var.HIend)
            self.var.YldFormCD = np.copy(self.var.YldForm)            
            self.var.FloweringEndCD = np.copy(self.var.FloweringEnd)  # only used in this function
            self.var.FloweringCD = np.copy(self.var.Flowering)

        # Pre-compute cumulative GDD during growing season
        if (self.var.CalendarType == 1 & self.var.SwitchGDD) | (self.var.CalendarType == 2):

            pd = np.copy(self.var.PlantingDate)
            hd = np.copy(self.var.HarvestDate)
            sd = self.var._modelTime.startTime.timetuple().tm_yday

            # NEW:
            isLeapYear1 = calendar.isleap(self.var._modelTime.startTime.year)
            isLeapYear2 = calendar.isleap(self.var._modelTime.startTime.year + 1)

            if isLeapYear1:
                pd[pd >= 60] += 1
                hd[(hd > pd) & (hd >= 60)] += 1

            if isLeapYear2:
                hd[(hd < pd) & (hd >= 60)] += 1
                
            hd[hd < pd] += 365  # already accounted for leap years so this should be OK

            cond = sd > pd
            pd[cond] += 365
            hd[cond] += 365

            # OLD:
            # # if start day of simulation is greater than planting day the
            # # first complete growing season will not be until the
            # # following year
            # pd[sd > pd] += 365
            # hd[sd > pd] += 365

            # # if start day is less than or equal to planting day, but
            # # harvest day is less than planting day, the harvest day will
            # # occur in the following year
            # hd[((sd <= pd) & (hd < pd))] += 365

            # # adjust values for leap year
            # isLeapYear1 = calendar.isleap(self.var._modelTime.startTime.year)
            # isLeapYear2 = calendar.isleap(self.var._modelTime.startTime.year + 1)
            # pd[(isLeapYear1 & (pd >= 60))] += 1  # TODO: check these
            # hd[(isLeapYear1 & (hd >= 60))] += 1
            # pd[(isLeapYear2 & (pd >= 425))] += 1
            # hd[(isLeapYear2 & (hd >= 425))] += 1            

            max_harvest_date = int(np.max(hd))
            day_idx = np.arange(sd, max_harvest_date + 1)[:,None,None,None] * np.ones_like(self.var.PlantingDate)[None,:,:,:]
            growing_season_idx = ((day_idx >= pd) & (day_idx <= hd))

            # Extract weather data for first growing season
            tmin = vos.netcdf2NumPyTimeSlice(self.var.tmpFileNC, self.var.tmnVarName,
                                             self.var._modelTime.startTime,
                                             self.var._modelTime.startTime + datetime.timedelta(int(max_harvest_date - sd)),
                                             cloneMapFileName = self.var.cloneMap,
                                             LatitudeLongitude = True)

            tmax = vos.netcdf2NumPyTimeSlice(self.var.tmpFileNC, self.var.tmxVarName,
                                             self.var._modelTime.startTime,
                                             self.var._modelTime.startTime + datetime.timedelta(int(max_harvest_date - sd)),
                                             cloneMapFileName = self.var.cloneMap,
                                             LatitudeLongitude = True)

            # broadcast to crop dimension
            tmax = tmax[:,None,:,:] * np.ones((self.var.PlantingDate.shape[0]))[None,:,None,None]
            tmin = tmin[:,None,:,:] * np.ones((self.var.PlantingDate.shape[0]))[None,:,None,None]

            # for convenience
            tupp = self.var.Tupp[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]
            tbase = self.var.Tbase[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]

            # calculate GDD according to the various methods
            if self.var.GDDmethod == 1:
                tmean = ((tmax + tmin) / 2)
                tmean = np.clip(tmean, self.var.Tbase, self.var.Tupp)
            elif self.var.GDDmethod == 2:
                tmax = np.clip(tmax, self.var.Tbase, self.var.Tupp)
                tmin = np.clip(tmin, self.var.Tbase, self.var.Tupp)
                tmean = ((tmax + tmin) / 2)
            elif self.var.GDDmethod == 3:
                tmax = np.clip(tmax, self.var.Tbase, self.var.Tupp)
                tmin = np.clip(tmin, None, self.var.Tupp)
                tmean = ((tmax + tmin) / 2)
                tmean = np.clip(tmean, self.var.Tbase, None)

            tmean[np.logical_not(growing_season_idx)] = 0
            tbase[np.logical_not(growing_season_idx)] = 0
            GDD = (tmean - tbase)
            GDDcum = np.cumsum(GDD, axis=0)

            # "Check if converting crop calendar to GDD mode"
            # if Mode == 1 & self.var.SwitchGDD:
            if self.var.CalendarType == 1 & self.var.SwitchGDD:

                # Find GDD equivalent for each crop calendar variable
                m,n,p = pd.shape  # crop,lat,lon
                I,J,K = np.ogrid[:m,:n,:p]

                emergence_idx = pd + self.var.EmergenceCD  # crop,lat,lon
                self.var.Emergence = GDDcum[emergence_idx,I,J,K]
                canopy10pct_idx = pd + self.var.Canopy10PctCD
                self.var.Canopy10Pct = GDDcum[canopy10pct_idx,I,J,K]
                maxrooting_idx = pd + self.var.MaxRootingCD
                self.var.MaxRooting = GDDcum[maxrooting_idx,I,J,K]
                maxcanopy_idx = pd + self.var.MaxCanopyCD
                self.var.MaxCanopy = GDDcum[maxcanopy_idx,I,J,K]
                canopydevend_idx = pd + self.var.CanopyDevEndCD
                self.var.CanopyDevEnd = GDDcum[canopydevend_idx,I,J,K]
                senescence_idx = pd + self.var.SenescenceCD
                self.var.Senescence = GDDcum[senescence_idx,I,J,K]
                maturity_idx = pd + self.var.MaturityCD
                self.var.Maturity = GDDcum[maturity_idx,I,J,K]
                histart_idx = pd + self.var.HIstartCD
                self.var.HIstart = GDDcum[histart_idx,I,J,K]
                hiend_idx = pd + self.var.HIendCD
                self.var.HIend = GDDcum[hiend_idx,I,J,K]
                yldform_idx = pd + self.var.YldFormCD
                self.var.YldForm = GDDcum[yldform_idx,I,J,K]

                cond2 = (self.var.CropType == 3)
                floweringend_idx = pd + self.var.FloweringEndCD
                self.var.FloweringEnd[cond2] = GDDcum[floweringend_idx,I,J,K][cond2]
                self.var.Flowering[cond2] = (self.var.FloweringEnd - self.var.HIstart)[cond2]

                # "Convert CGC to GDD mode"
                # self.var.CGC_CD = self.var.CGC
                self.var.CGC = (np.log((((0.98 * self.var.CCx) - self.var.CCx) * self.var.CC0) / (-0.25 * (self.var.CCx ** 2)))) / (-(self.var.MaxCanopy - self.var.Emergence))

                # "Convert CDC to GDD mode"
                # self.var.CDC_CD = self.var.CDC
                tCD = self.var.MaturityCD - self.var.SenescenceCD
                tCD[tCD <= 0] = 1
                tGDD = self.var.Maturity - self.var.Senescence
                tGDD[tGDD <= 0] = 5
                self.var.CDC = (self.var.CCx / tGDD) * np.log(1 + ((1 - self.var.CCi / self.var.CCx) / 0.05))

                # "Set calendar type to GDD mode"
                self.var._configuration.cropOptions['CalendarType'] = "2"
            
            # elif Mode == 2:
            elif self.var.CalendarType == 2:

                # "Find calendar days [equivalent] for some variables"

                # "1 Calendar days from sowing to maximum canopy cover"

                # # TODO: check this indexing
                # day_idx = np.arange(0, GDD.shape[0])[:,None,None,None] * np.ones((self.var.nCrop, self.var.nLon, self.var.nLat))[None,:,:,:]
                # pd,hd = self.var.adjust_planting_and_harvesting_date(self.var.modelTime.startTime)

                maxcanopy_idx = np.copy(day_idx)
                maxcanopy_idx[np.logical_not(GDDcum > self.var.MaxCanopy)] = 999
                # maxcanopy_idx[np.logical_not(GDDcum > self.var.MaxCanopy)] = np.nan
                maxcanopy_idx = np.nanmin(maxcanopy_idx, axis=0)
                self.var.MaxCanopyCD = maxcanopy_idx - pd + 1

                # "2 Calendar days from sowing to end of vegetative growth"
                canopydevend_idx = np.copy(day_idx)
                canopydevend_idx[np.logical_not(GDDcum > self.var.CanopyDevEnd)] = 999
                # canopydevend_idx[np.logical_not(GDDcum > self.var.CanopyDevEnd)] = np.nan
                canopydevend_idx = np.nanmin(canopydevend_idx, axis=0)
                self.var.CanopyDevEndCD = canopydevend_idx - pd + 1

                # "3 Calendar days from sowing to start of yield formation"
                histart_idx = np.copy(day_idx)
                histart_idx[np.logical_not(GDDcum > self.var.HIstart)] = 999
                # histart_idx[np.logical_not(GDDcum > self.var.HIstart)] = np.nan
                histart_idx = np.nanmin(histart_idx, axis=0)
                self.var.HIstartCD = histart_idx - pd + 1

                # "4 Calendar days from sowing to end of yield formation"
                hiend_idx = np.copy(day_idx)
                hiend_idx[np.logical_not(GDDcum > self.var.HIend)] = 999
                # hiend_idx[np.logical_not(GDDcum > self.var.HIend)] = np.nan
                hiend_idx = np.nanmin(hiend_idx, axis=0)
                self.var.HIendCD = hiend_idx - pd + 1

                # "Duration of yield formation in calendar days"
                self.var.YldFormCD = self.var.HIendCD - self.var.HIstartCD

                cond1 = (self.var.CropType == 3)

                # "1 Calendar days from sowing to end of flowering"
                floweringend_idx = np.copy(day_idx)
                floweringend_idx[np.logical_not(GDDcum > self.var.FloweringEnd)] = 999
                # floweringend_idx[np.logical_not(GDDcum > self.var.FloweringEnd)] = np.nan
                floweringend_idx = np.nanmin(floweringend_idx, axis=0)
                FloweringEnd = floweringend_idx - pd + 1

                # "2 Duration of flowering in calendar days"
                self.var.FloweringCD[cond1] = (FloweringEnd - self.var.HIstartCD)[cond1]
                
    def update_crop_parameters(self):
        """Function to update certain crop parameters for current 
        time step (equivalent to lines 97-163 in 
        compute_crop_calendar)
        """
        pd = np.copy(self.var.PlantingDateAdj)
        hd = np.copy(self.var.HarvestDateAdj)
        hd[hd < pd] += 365
        sd = self.var._modelTime.currTime.timetuple().tm_yday

        # Update certain crop parameters if using GDD mode
        if (self.var.CalendarType == 2):

            cond1 = ((self.var.GrowingSeason) & (self.var._modelTime.doy == pd))
            pd[np.logical_not(cond1)] = 0
            hd[np.logical_not(cond1)] = 0
            max_harvest_date = int(np.max(hd))

            if (max_harvest_date > 0):

                # Dimension (day,crop,lat,lon)
                day_idx = np.arange(sd, max_harvest_date + 1)[:,None,None,None] * np.ones_like(self.var.PlantingDate)[None,:,:,:]
                growing_season_idx = ((day_idx >= pd) & (day_idx <= hd))

                # Extract weather data for first growing season
                tmin = vos.netcdf2NumPyTimeSlice(self.var.tmpFileNC, self.var.tmnVarName,
                                                 self.var._modelTime.currTime,
                                                 self.var._modelTime.currTime + datetime.timedelta(int(max_harvest_date - sd)),
                                                 cloneMapFileName = self.var.cloneMap,
                                                 LatitudeLongitude = True)

                tmax = vos.netcdf2NumPyTimeSlice(self.var.tmpFileNC, self.var.tmxVarName,
                                                 self.var._modelTime.currTime,
                                                 self.var._modelTime.currTime + datetime.timedelta(int(max_harvest_date - sd)),
                                                 cloneMapFileName = self.var.cloneMap,
                                                 LatitudeLongitude = True)

                # broadcast to crop dimension
                tmax = tmax[:,None,:,:] * np.ones((self.var.nCrop))[None,:,None,None]
                tmin = tmin[:,None,:,:] * np.ones((self.var.nCrop))[None,:,None,None]

                # for convenience
                tupp = self.var.Tupp[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]
                tbase = self.var.Tbase[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]

                # calculate GDD according to the various methods
                if self.var.GDDmethod == 1:
                    tmean = ((tmax + tmin) / 2)
                    tmean = np.clip(tmean, self.var.Tbase, self.var.Tupp)
                elif self.var.GDDmethod == 2:
                    tmax = np.clip(tmax, self.var.Tbase, self.var.Tupp)
                    tmin = np.clip(tmin, self.var.Tbase, self.var.Tupp)
                    tmean = ((tmax + tmin) / 2)
                elif self.var.GDDmethod == 3:
                    tmax = np.clip(tmax, self.var.Tbase, self.var.Tupp)
                    tmin = np.clip(tmin, None, self.var.Tupp)
                    tmean = ((tmax + tmin) / 2)
                    tmean = np.clip(tmean, self.var.Tbase, None)

                tmean[np.logical_not(growing_season_idx)] = 0
                tbase[np.logical_not(growing_season_idx)] = 0
                GDD = (tmean - tbase)
                GDDcum = np.cumsum(GDD, axis=0)
        
                # 1 - Calendar days from sowing to maximum canopy cover
                maxcanopy_idx = np.copy(day_idx)
                maxcanopy_idx[np.logical_not(GDDcum > self.var.MaxCanopy)] = 999
                # maxcanopy_idx[np.logical_not(GDDcum > self.var.MaxCanopy)] = np.nan
                maxcanopy_idx = np.nanmin(maxcanopy_idx, axis=0)
                MaxCanopyCD = (maxcanopy_idx - pd + 1)
                self.var.MaxCanopyCD[cond1] = MaxCanopyCD[cond1]

                # 2 - Calendar days from sowing to end of vegetative growth
                canopydevend_idx = np.copy(day_idx)
                canopydevend_idx[np.logical_not(GDDcum > self.var.CanopyDevEnd)] = 999
                # canopydevend_idx[np.logical_not(GDDcum > self.var.CanopyDevEnd)] = np.nan
                canopydevend_idx = np.nanmin(canopydevend_idx, axis=0)
                CanopyDevEndCD = canopydevend_idx - pd + 1
                self.var.CanopyDevEndCD[cond1] = CanopyDevEndCD[cond1]

                # 3 - Calendar days from sowing to start of yield formation
                histart_idx = np.copy(day_idx)
                histart_idx[np.logical_not(GDDcum > self.var.HIstart)] = 999
                # histart_idx[np.logical_not(GDDcum > self.var.HIstart)] = np.nan
                histart_idx = np.nanmin(histart_idx, axis=0)
                HIstartCD = histart_idx - pd + 1
                self.var.HIstartCD[cond1] = HIstartCD[cond1]

                # 4 - Calendar days from sowing to end of yield formation
                hiend_idx = np.copy(day_idx)
                hiend_idx[np.logical_not(GDDcum > self.var.HIend)] = 999
                # hiend_idx[np.logical_not(GDDcum > self.var.HIend)] = np.nan
                hiend_idx = np.nanmin(hiend_idx, axis=0)
                HIendCD = hiend_idx - pd + 1
                self.var.HIendCD[cond1] = HIendCD[cond1]

                # Duration of yield formation in calendar days
                self.var.YldFormCD[cond1] = (self.var.HIendCD - self.var.HIstartCD)[cond1]

                cond11 = (cond1 & (self.var.CropType == 3))

                # 1 Calendar days from sowing to end of flowering
                floweringend_idx = np.copy(day_idx)
                floweringend_idx[np.logical_not(GDDcum > self.var.FloweringEnd)] = 999
                # floweringend_idx[np.logical_not(GDDcum > self.var.FloweringEnd)] = np.nan
                floweringend_idx = np.nanmin(floweringend_idx, axis=0)
                FloweringEnd = floweringend_idx - pd + 1

                # 2 Duration of flowering in calendar days
                self.var.FloweringCD[cond11] = (FloweringEnd - self.var.HIstartCD)[cond11]

                # Harvest index growth coefficient
                self.calculate_HIGC()

                # Days to linear HI switch point
                self.calculate_HI_linear()

    def dynamic(self):
        """Function to update parameters for current crop grown as well 
        as counters pertaining to crop growth
        """
        # Update crop parameters for currently grown crops
        # self.compute_water_productivity_adjustment_factor()
        self.adjust_planting_and_harvesting_date()
        self.update_growing_season()
        self.compute_water_productivity_adjustment_factor()
        self.update_crop_parameters()
        
class FAO56CropParameters(CropParameters):

    def initial(self):

        super(FAO56CropParameters, self).initial()
        
        # Declare variables
        self.var.crop_parameters_to_read = [
            'PlantingDate','HarvestDate',
            'L_ini','L_dev','L_mid','L_late','Kc_ini','Kc_mid','Kc_end','p_std',
            'Zmin','Zmax','Ky']

        self.var.crop_parameters_to_compute = [
            'L_ini_day','L_dev_day','L_mid_day','L_late_day',
            'PlantingDateAdj','HarvestDateAdj']

        # initialise parameters
        self.var.crop_parameter_names = self.var.crop_parameters_to_read + self.var.crop_parameters_to_compute
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        for param in self.var.crop_parameter_names:
            vars(self.var)[param] = np.copy(arr_zeros)

        # potential yield
        self.var.crop_parameter_names += ['Yx']
        self.var.PotYieldFileNC = self.var._configuration.cropOptions['PotYieldNC']
        self.var.PotYieldVarName = 'Yx' 
        if 'PotYieldVariableName' in self.var._configuration.cropOptions:
            self.var.PotYieldVarName = self.var._configuration.cropOptions['PotYieldVariableName']
        # self.var.co2_set_per_year  = False

        # find out whether potential crop yield parameter is dynamic
        self.var.AnnualChangeInPotYield = False
        if 'AnnualChangeInPotYield' in self.var._configuration.cropOptions:
            if self.var.AnnualChangeInPotYield == "True":
                self.var.AnnualChangeInPotYield = True
                
        # # TODO: only read data if the year has changed        
        # date = '%04i-%02i-%02i' %(self.var._modelTime.year, 1, 1)
        # self.var.Yx = vos.netcdf2PCRobjClone(self.var.PotYieldFileNC,
        #                                      self.var.PotYieldVarName,
        #                                      date,
        #                                      useDoy = None,
        #                                      cloneMapFileName = self.var.cloneMap,
        #                                      LatitudeLongitude = True)
         
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        self.var.CropDead = np.copy(arr_zeros).astype(bool)
        self.var.CropMature = np.copy(arr_zeros).astype(bool)        
        self.read() 

    def compute_growth_stage_length(self):
        nday = self.var.HarvestDateAdj - self.var.PlantingDateAdj
        nday[nday < 0] += 365
        self.var.L_ini_day = np.round(self.var.L_ini * nday)
        self.var.L_dev_day = np.round(self.var.L_dev * nday)
        self.var.L_mid_day = np.round(self.var.L_mid * nday)
        self.var.L_late_day = np.round(self.var.L_late * nday)  # TODO

    def read_potential_crop_yield(self):
        if self.var.AnnualChangeInPotYield:
            if self.var._modelTime.timeStepPCR == 1 or self.var._modelTime.doy == 1:
                date = '%04i-%02i-%02i' %(self.var._modelTime.year, 1, 1)
                self.var.Yx = vos.netcdf2PCRobjClone(
                    self.var.PotYieldFileNC,
                    self.var.PotYieldVarName,
                    date,
                    useDoy = None,
                    cloneMapFileName = self.var.cloneMap,
                    LatitudeLongitude = True)
        else:
            if self.var._modelTime.timeStepPCR == 1:
                self.var.Yx = vos.netcdf2PCRobjCloneWithoutTime(
                    self.var.PotYieldFileNC,
                    self.var.PotYieldVarName,
                    cloneMapFileName = self.var.cloneMap)

    def dynamic(self):
        self.adjust_planting_and_harvesting_date()
        self.compute_growth_stage_length()
        self.update_growing_season()
        self.read_potential_crop_yield()
