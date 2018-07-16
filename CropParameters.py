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
        self.var.nRotation = int(self.var._configuration.cropOptions['nRotation'])
        self.var.cropParameterFileNC = str(self.var._configuration.cropOptions['cropParameterNC'])
        self.var.crop_parameters_to_read = []
        self.var.crop_parameters_to_compute = []

    def read(self):
        """Function to read crop input parameters"""
        if len(self.var.crop_parameters_to_read) > 0:
            for param in self.var.crop_parameters_to_read:
                nm = '_' + param
                vars(self.var)[nm] = vos.netcdf2PCRobjCloneWithoutTime(
                    self.var.cropParameterFileNC,
                    param,
                    cloneMapFileName=self.var.cloneMap)

    def read_crop_sequence(self):
        """Function to read crop sequence information"""        
        # Matrix showing specific crop rotations
        CropSequence = vos.netcdf2PCRobjCloneWithoutTime(
            self.var.cropParameterFileNC,
            'CropSequence',
            cloneMapFileName=self.var.cloneMap)

        CropSequence = (
            CropSequence[:,:,None,None]
            * np.ones((self.var.nLat, self.var.nLon))[None,None,:,:])
        self.var.CropSequence = CropSequence.astype(bool)

    def adjust_planting_and_harvesting_date(self):
        pd = np.copy(self.var._PlantingDate)
        hd = np.copy(self.var._HarvestDate)
        st = self.var._modelTime.currTime
        sd = self.var._modelTime.currTime.timetuple().tm_yday

        # adjust values for leap year (objective is to preserve date)
        isLeapYear1 = calendar.isleap(st.year)
        isLeapYear2 = calendar.isleap(st.year + 1)
        pd[(isLeapYear1 & (pd >= 60))] += 1  # TODO: check these
        hd[(isLeapYear1 & (hd >= 60) & (hd > pd))] += 1
        hd[(isLeapYear2 & (hd >= 60) & (hd < pd))] += 1

        self.var._PlantingDate = np.copy(pd)
        self.var._HarvestDate = np.copy(hd)
        
        # if harvest day is less than planting day, harvest will take place
        # in following year (leap year already accounted for, so add 365)
        hd[hd < pd] += 365
        self.var._PlantingDateAdj = pd
        self.var._HarvestDateAdj = hd

    def update_growing_season(self):
        cond1 = ((self.var._PlantingDateAdj <= self.var._modelTime.doy) & (self.var._modelTime.doy <= self.var._HarvestDateAdj))
        hd_arr = (np.datetime64(str(datetime.datetime(self.var._modelTime.year, 1, 1))) + np.array(self.var._HarvestDateAdj - 1, dtype='timedelta64[D]'))
        cond3 = hd_arr <= np.datetime64(str(self.var._modelTime.endTime))        
        self.var.GrowingSeason = ((cond1 & cond3))
        # TODO: when rotation dimension is removed, compute GrowingSeasonDayOne here, and delete GrowingSeason class
        
    def select_crop_parameters(self):
        """Function to select parameters for currently grown crops"""
        
        # Add a rotation dimension to GrowingSeason array and multiply it
        # by CropSequence, with the resulting array showing which crop (if any)
        # is currently grown in each rotation considered.
        GrowingSeason = self.var.GrowingSeason[:,None,:,:] * np.ones((self.var.nRotation))[None,:,None,None]
        GrowingSeason *= self.var.CropSequence  # crop,rotation,lat,lon
        
        # Lastly, check whether the crop has died or reached maturity, then
        # see which rotations have crops currently growing
        self.var.GrowingSeasonIndex = np.any(GrowingSeason, axis=0)
        
        # Get index of crops currently grown
        CropIndex = (np.arange(0, self.var.nCrop)[:,None,None,None] * np.ones((self.var.nRotation, self.var.nLon, self.var.nLat))[None,:,:,:])
        CropIndex *= GrowingSeason
        self.var.CropIndex = np.max(CropIndex, axis=0).astype(int)
        
        # A CropIndex value of 0 currently means either that the crop is not
        # currently grown, or that the first crop (corresponding to index 0)
        # is grown. This confusion is handled in the next section by multiplying
        # the parameter value by GrowingSeason, such that it has a value of
        # zero if no crops are growing.

        # Select crop parameters for current day
        I,J,K = np.ogrid[:self.var.nRotation,:self.var.nLat,:self.var.nLon]
        for nm in self.var.crop_parameter_names:
            param = getattr(self.var, '_' + nm)
            param = param[:,None,:,:] * np.ones((self.var.nRotation))[None,:,None,None]
            vars(self.var)[nm] = param[self.var.CropIndex,I,J,K] * self.var.GrowingSeasonIndex
        
class AQCropParameters(CropParameters):

    def initial(self):

        # TODO: move these to main model class?
        # self.var.nCrop = int(self.var._configuration.cropOptions['nCrop'])
        # self.var.nRotation = int(self.var._configuration.cropOptions['nRotation'])
        # self.var.cropParameterFileNC = str(self.var._configuration.cropOptions['cropParameterNC'])
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
            'PlantingDateAdj','HarvestDateAdj']  # TODO: CGC and CDC are in both list - try removing them here?

        self.var.crop_parameter_names = self.var.crop_parameters_to_read + self.var.crop_parameters_to_compute
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        for param in self.var.crop_parameter_names:
            nm = '_' + param
            vars(self.var)[nm] = arr_zeros
            
        # self.var.CropSequence = np.zeros((self.var.nCrop, self.var.nRotation, self.var.nLat, self.var.nLon))
        self.read()
        self.read_crop_sequence()
        self.compute_crop_parameters()

    def compute_crop_parameters(self):
        """Function to compute additional crop variables required to run 
        AquaCrop"""
        
        # The following adapted from lines 160-240 of AOS_ComputeVariables.m

        # Fractional canopy cover size at emergence
        self.var._CC0 = np.round(10000 * self.var._PlantPop * self.var._SeedSize * 10 ** -8) / 10000
        
        # Root extraction terms
        S1 = np.copy(self.var._SxTopQ)
        S2 = np.copy(self.var._SxBotQ)
        self.var._SxTop = np.zeros_like(self.var._SxTopQ)
        self.var._SxBot = np.zeros_like(self.var._SxBotQ)
        # self.var.SxTop = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        # self.var.SxBot = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))

        cond1 = (S1 == S2)
        self.var._SxTop[cond1] = S1[cond1]
        self.var._SxBot[cond1] = S2[cond1]

        cond2 = np.logical_not(cond1)
        cond21 = (cond2 & (self.var._SxTopQ < self.var._SxBotQ))
        S1[cond21] = self.var._SxBotQ[cond21]
        S2[cond21] = self.var._SxTopQ[cond21]
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

        cond24 = (cond2 & (self.var._SxTopQ > self.var._SxBotQ))
        self.var._SxTop[cond24] = SS1[cond24]
        self.var._SxBot[cond24] = SS2[cond24]

        cond25 = (cond2 & np.logical_not(cond24))
        self.var._SxTop[cond25] = SS2[cond25]
        self.var._SxBot[cond25] = SS1[cond25]

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
                                           * self.var._bsted + fw
                                           * ((self.var._bsted * self.var._fsink)
                                              + (self.var._bface
                                                 * (1 - self.var._fsink))))))

        # Consider crop type
        ftype = (40 - self.var._WP) / (40 - 20)
        ftype = np.clip(ftype, 0, 1)
        self.var._fCO2 = 1 + ftype * (fCO2 - 1)
        
    def calculate_HI_linear(self):
        """Function to calculate time to switch to linear harvest index 
        build-up, and associated linear rate of build-up. Only for 
        fruit/grain crops
        """
        # Determine linear switch point
        cond1 = (self.var._CropType == 3)
        ti = np.zeros_like(self.var._CropType)
        # ti = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        tmax = self.var._YldFormCD
        HIest = np.zeros_like(self.var._HIini)
        # HIest = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        HIprev = self.var._HIini

        # Iterate to find linear switch point
        while np.any((self.var._CropType == 3) & (HIest <= self.var._HI0) & (ti < tmax)):    
            ti += 1
            HInew = ((self.var._HIini * self.var._HI0) / (self.var._HIini + (self.var._HI0 - self.var._HIini) * np.exp(-self.var._HIGC * ti)))
            HIest = (HInew + (tmax - ti) * (HInew - HIprev))
            HIprev = HInew

        self.var._tLinSwitch = ti
        # self.var.tLinSwitch[self.var.CropType != 3] = np.nan
            
        # Determine linear build-up rate
        cond1 = (self.var._tLinSwitch > 0)
        HIest[cond1] = ((self.var._HIini * self.var._HI0) / (self.var._HIini + (self.var._HI0 - self.var._HIini) * np.exp(-self.var._HIGC * self.var._tLinSwitch)))[cond1]
        HIest[np.logical_not(cond1)] = 0
        self.var._dHILinear = ((self.var._HI0 - HIest) / (tmax - self.var._tLinSwitch))  # dHILin will be set to nan in the same cells as tSwitch
                    
    def calculate_HIGC(self):
        """Function to calculate harvest index growth coefficient"""
        # Total yield formation days
        tHI = np.copy(self.var._YldFormCD)

        # Iteratively estimate HIGC
        self.var._HIGC = np.full((self.var._HIini.shape), 0.001)
        # self.var.HIGC = np.full((self.var.nCrop, self.var.nLat, self.var.nLon), 0.001)
        HIest = np.zeros_like(self.var._HIini)
        # HIest = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        while np.any(HIest < (0.98 * self.var._HI0)):
            cond1 = (HIest < (0.98 * self.var._HI0))
            self.var._HIGC += 0.001
            HIest = ((self.var._HIini * self.var._HI0) / (self.var._HIini + (self.var._HI0 - self.var._HIini) * np.exp(-self.var._HIGC * tHI)))
            
        self.var._HIGC[HIest >= self.var._HI0] -= 0.001

    def compute_crop_calendar(self):
       
        # "Time from sowing to end of vegetative growth period"
        cond1 = (self.var._Determinant == 1)
        self.var._CanopyDevEnd = np.copy(self.var._Senescence)
        self.var._CanopyDevEnd[cond1] = (np.round(self.var._HIstart + (self.var._Flowering / 2)))[cond1]

        # "Time from sowing to 10% canopy cover (non-stressed conditions)
        self.var._Canopy10Pct = np.round(self.var._Emergence + (np.log(0.1 / self.var._CC0) / self.var._CGC))

        # "Time from sowing to maximum canopy cover (non-stressed conditions)
        self.var._MaxCanopy = np.round(self.var._Emergence + (np.log((0.25 * self.var._CCx * self.var._CCx / self.var._CC0) / (self.var._CCx - (0.98 * self.var._CCx))) / self.var._CGC))

        # "Time from sowing to end of yield formation"
        self.var._HIend = self.var._HIstart + self.var._YldForm

        cond2 = (self.var._CropType == 3)

        # TODO: declare these in __init__
        arr_zeros = np.zeros_like(self.var._CropType)
        self.var._FloweringEnd = arr_zeros 
        self.var._FloweringEndCD = arr_zeros
        self.var._FloweringCD = arr_zeros
                
        self.var._FloweringEnd[cond2] = (self.var._HIstart + self.var._Flowering)[cond2]
        self.var._FloweringEndCD[cond2] = self.var._FloweringEndCD[cond2]
        self.var._FloweringCD[cond2] = self.var._Flowering[cond2]
        
        # Mode = self.var.CalendarType
        # if Mode == 1:
        if self.var.CalendarType == 1:
            # "Duplicate calendar values (needed to minimise if-statements when
            # switching between GDD and CD runs)
            self.var._EmergenceCD = np.copy(self.var._Emergence)      # only used in this function
            self.var._Canopy10PctCD = np.copy(self.var._Canopy10Pct)  # only used in this function
            self.var._MaxRootingCD = np.copy(self.var._MaxRooting)    # only used in this function
            self.var._SenescenceCD = np.copy(self.var._Senescence)    # only used in this function
            self.var._MaturityCD = np.copy(self.var._Maturity)        # only used in this function
            self.var._MaxCanopyCD = np.copy(self.var._MaxCanopy)
            self.var._CanopyDevEndCD = np.copy(self.var._CanopyDevEnd)
            self.var._HIstartCD = np.copy(self.var._HIstart)
            self.var._HIendCD = np.copy(self.var._HIend)
            self.var._YldFormCD = np.copy(self.var._YldForm)            
            self.var._FloweringEndCD = np.copy(self.var._FloweringEnd)  # only used in this function
            self.var._FloweringCD = np.copy(self.var._Flowering)

        # Pre-compute cumulative GDD during growing season
        if (self.var.CalendarType == 1 & self.var.SwitchGDD) | (self.var.CalendarType == 2):

            pd = np.copy(self.var._PlantingDate)
            hd = np.copy(self.var._HarvestDate)
            sd = self.var._modelTime.startTime.timetuple().tm_yday

            # if start day of simulation is greater than planting day the
            # first complete growing season will not be until the
            # following year
            pd[sd > pd] += 365
            hd[sd > pd] += 365

            # if start day is less than or equal to planting day, but
            # harvest day is less than planting day, the harvest day will
            # occur in the following year
            hd[((sd <= pd) & (hd < pd))] += 365

            # adjust values for leap year
            isLeapYear1 = calendar.isleap(self.var._modelTime.startTime.year)
            isLeapYear2 = calendar.isleap(self.var._modelTime.startTime.year + 1)
            pd[(isLeapYear1 & (pd >= 60))] += 1  # TODO: check these
            hd[(isLeapYear1 & (hd >= 60))] += 1
            pd[(isLeapYear2 & (pd >= 425))] += 1
            hd[(isLeapYear2 & (hd >= 425))] += 1            

            max_harvest_date = int(np.max(hd))
            day_idx = np.arange(sd, max_harvest_date + 1)[:,None,None,None] * np.ones_like(self.var._PlantingDate)[None,:,:,:]
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
            tmax = tmax[:,None,:,:] * np.ones((self.var._PlantingDate.shape[0]))[None,:,None,None]
            tmin = tmin[:,None,:,:] * np.ones((self.var._PlantingDate.shape[0]))[None,:,None,None]
            # tmax = tmax[:,None,:,:] * np.ones((self.var.nCrop))[None,:,None,None]
            # tmin = tmin[:,None,:,:] * np.ones((self.var.nCrop))[None,:,None,None]

            # for convenience
            tupp = self.var._Tupp[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]
            tbase = self.var._Tbase[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]

            # calculate GDD according to the various methods
            if self.var.GDDmethod == 1:
                tmean = ((tmax + tmin) / 2)
                tmean = np.clip(tmean, self.var._Tbase, self.var._Tupp)
            elif self.var.GDDmethod == 2:
                tmax = np.clip(tmax, self.var._Tbase, self.var._Tupp)
                tmin = np.clip(tmin, self.var._Tbase, self.var._Tupp)
                tmean = ((tmax + tmin) / 2)
            elif self.var.GDDmethod == 3:
                tmax = np.clip(tmax, self.var._Tbase, self.var._Tupp)
                tmin = np.clip(tmin, None, self.var._Tupp)
                tmean = ((tmax + tmin) / 2)
                tmean = np.clip(tmean, self.var._Tbase, None)

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

                emergence_idx = pd + self.var._EmergenceCD  # crop,lat,lon
                self.var._Emergence = GDDcum[emergence_idx,I,J,K]
                canopy10pct_idx = pd + self.var._Canopy10PctCD
                self.var._Canopy10Pct = GDDcum[canopy10pct_idx,I,J,K]
                maxrooting_idx = pd + self.var._MaxRootingCD
                self.var._MaxRooting = GDDcum[maxrooting_idx,I,J,K]
                maxcanopy_idx = pd + self.var._MaxCanopyCD
                self.var._MaxCanopy = GDDcum[maxcanopy_idx,I,J,K]
                canopydevend_idx = pd + self.var._CanopyDevEndCD
                self.var._CanopyDevEnd = GDDcum[canopydevend_idx,I,J,K]
                senescence_idx = pd + self.var._SenescenceCD
                self.var._Senescence = GDDcum[senescence_idx,I,J,K]
                maturity_idx = pd + self.var._MaturityCD
                self.var._Maturity = GDDcum[maturity_idx,I,J,K]
                histart_idx = pd + self.var._HIstartCD
                self.var._HIstart = GDDcum[histart_idx,I,J,K]
                hiend_idx = pd + self.var._HIendCD
                self.var._HIend = GDDcum[hiend_idx,I,J,K]
                yldform_idx = pd + self.var._YldFormCD
                self.var._YldForm = GDDcum[yldform_idx,I,J,K]

                cond2 = (self.var._CropType == 3)
                floweringend_idx = pd + self.var._FloweringEndCD
                self.var._FloweringEnd[cond2] = GDDcum[floweringend_idx,I,J,K][cond2]
                self.var._Flowering[cond2] = (self.var._FloweringEnd - self.var._HIstart)[cond2]

                # "Convert CGC to GDD mode"
                # self.var.CGC_CD = self.var.CGC
                self.var._CGC = (np.log((((0.98 * self.var._CCx) - self.var._CCx) * self.var._CC0) / (-0.25 * (self.var._CCx ** 2)))) / (-(self.var._MaxCanopy - self.var._Emergence))

                # "Convert CDC to GDD mode"
                # self.var.CDC_CD = self.var.CDC
                tCD = self.var._MaturityCD - self.var._SenescenceCD
                tCD[tCD <= 0] = 1
                tGDD = self.var._Maturity - self.var._Senescence
                tGDD[tGDD <= 0] = 5
                self.var._CDC = (self.var._CCx / tGDD) * np.log(1 + ((1 - self.var._CCi / self.var._CCx) / 0.05))

                # "Set calendar type to GDD mode"
                self.var._configuration.cropOptions['CalendarType'] = "2"
            
            # elif Mode == 2:
            elif self.var.CalendarType == 2:

                # "Find calendar days [equivalent] for some variables"

                # "1 Calendar days from sowing to maximum canopy cover"

                # # TODO: check this indexing
                # day_idx = np.arange(0, GDD.shape[0])[:,None,None,None] * np.ones((self.var.nCrop, self.var.nLon, self.var.nLat))[None,:,:,:]
                # pd,hd = self.var.adjust_planting_and_harvesting_date(self.var._modelTime.startTime)

                maxcanopy_idx = np.copy(day_idx)
                maxcanopy_idx[np.logical_not(GDDcum > self.var._MaxCanopy)] = 999
                # maxcanopy_idx[np.logical_not(GDDcum > self.var.MaxCanopy)] = np.nan
                maxcanopy_idx = np.nanmin(maxcanopy_idx, axis=0)
                self.var._MaxCanopyCD = maxcanopy_idx - pd + 1

                # "2 Calendar days from sowing to end of vegetative growth"
                canopydevend_idx = np.copy(day_idx)
                canopydevend_idx[np.logical_not(GDDcum > self.var._CanopyDevEnd)] = 999
                # canopydevend_idx[np.logical_not(GDDcum > self.var.CanopyDevEnd)] = np.nan
                canopydevend_idx = np.nanmin(canopydevend_idx, axis=0)
                self.var._CanopyDevEndCD = canopydevend_idx - pd + 1

                # "3 Calendar days from sowing to start of yield formation"
                histart_idx = np.copy(day_idx)
                histart_idx[np.logical_not(GDDcum > self.var._HIstart)] = 999
                # histart_idx[np.logical_not(GDDcum > self.var.HIstart)] = np.nan
                histart_idx = np.nanmin(histart_idx, axis=0)
                self.var._HIstartCD = histart_idx - pd + 1

                # "4 Calendar days from sowing to end of yield formation"
                hiend_idx = np.copy(day_idx)
                hiend_idx[np.logical_not(GDDcum > self.var._HIend)] = 999
                # hiend_idx[np.logical_not(GDDcum > self.var.HIend)] = np.nan
                hiend_idx = np.nanmin(hiend_idx, axis=0)
                self.var._HIendCD = hiend_idx - pd + 1

                # "Duration of yield formation in calendar days"
                self.var._YldFormCD = self.var._HIendCD - self.var._HIstartCD

                cond1 = (self.var._CropType == 3)

                # "1 Calendar days from sowing to end of flowering"
                floweringend_idx = np.copy(day_idx)
                floweringend_idx[np.logical_not(GDDcum > self.var._FloweringEnd)] = 999
                # floweringend_idx[np.logical_not(GDDcum > self.var.FloweringEnd)] = np.nan
                floweringend_idx = np.nanmin(floweringend_idx, axis=0)
                FloweringEnd = floweringend_idx - pd + 1

                # "2 Duration of flowering in calendar days"
                self.var._FloweringCD[cond1] = (FloweringEnd - self.var._HIstartCD)[cond1]

    def update_crop_parameters(self):
        """Function to update certain crop parameters for current 
        time step (equivalent to lines 97-163 in 
        compute_crop_calendar)
        """
        pd = np.copy(self.var._PlantingDateAdj)
        hd = np.copy(self.var._HarvestDateAdj)
        sd = self.var._modelTime.currTime.timetuple().tm_yday

        # Update certain crop parameters if using GDD mode
        if (self.var.CalendarType == 2):

            cond1 = ((self.var.GrowingSeason) & (self.var._modelTime.doy == pd))
            pd[np.logical_not(cond1)] = 0
            hd[np.logical_not(cond1)] = 0
            max_harvest_date = int(np.max(hd))

            if (max_harvest_date > 0):
                
                day_idx = np.arange(sd, max_harvest_date + 1)[:,None,None,None] * np.ones_like(self.var._PlantingDate)[None,:,:,:]
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
                tupp = self.var._Tupp[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]
                tbase = self.var._Tbase[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]

                # calculate GDD according to the various methods
                if self.var.GDDmethod == 1:
                    tmean = ((tmax + tmin) / 2)
                    tmean = np.clip(tmean, self.var._Tbase, self.var._Tupp)
                elif self.var.GDDmethod == 2:
                    tmax = np.clip(tmax, self.var._Tbase, self.var._Tupp)
                    tmin = np.clip(tmin, self.var._Tbase, self.var._Tupp)
                    tmean = ((tmax + tmin) / 2)
                elif self.var.GDDmethod == 3:
                    tmax = np.clip(tmax, self.var._Tbase, self.var._Tupp)
                    tmin = np.clip(tmin, None, self.var._Tupp)
                    tmean = ((tmax + tmin) / 2)
                    tmean = np.clip(tmean, self.var._Tbase, None)

                tmean[np.logical_not(growing_season_idx)] = 0
                tbase[np.logical_not(growing_season_idx)] = 0
                GDD = (tmean - tbase)
                GDDcum = np.cumsum(GDD, axis=0)
        
                # 1 - Calendar days from sowing to maximum canopy cover
                maxcanopy_idx = np.copy(day_idx)
                maxcanopy_idx[np.logical_not(GDDcum > self.var._MaxCanopy)] = 999
                # maxcanopy_idx[np.logical_not(GDDcum > self.var._MaxCanopy)] = np.nan
                maxcanopy_idx = np.nanmin(maxcanopy_idx, axis=0)
                MaxCanopyCD = (maxcanopy_idx - pd + 1)
                self.var._MaxCanopyCD[cond1] = MaxCanopyCD[cond1]

                # 2 - Calendar days from sowing to end of vegetative growth
                canopydevend_idx = np.copy(day_idx)
                canopydevend_idx[np.logical_not(GDDcum > self.var._CanopyDevEnd)] = 999
                # canopydevend_idx[np.logical_not(GDDcum > self.var._CanopyDevEnd)] = np.nan
                canopydevend_idx = np.nanmin(canopydevend_idx, axis=0)
                CanopyDevEndCD = canopydevend_idx - pd + 1
                self.var._CanopyDevEndCD[cond1] = CanopyDevEndCD[cond1]

                # 3 - Calendar days from sowing to start of yield formation
                histart_idx = np.copy(day_idx)
                histart_idx[np.logical_not(GDDcum > self.var._HIstart)] = 999
                # histart_idx[np.logical_not(GDDcum > self.var._HIstart)] = np.nan
                histart_idx = np.nanmin(histart_idx, axis=0)
                HIstartCD = histart_idx - pd + 1
                self.var._HIstartCD[cond1] = HIstartCD[cond1]

                # 4 - Calendar days from sowing to end of yield formation
                hiend_idx = np.copy(day_idx)
                hiend_idx[np.logical_not(GDDcum > self.var._HIend)] = 999
                # hiend_idx[np.logical_not(GDDcum > self.var._HIend)] = np.nan
                hiend_idx = np.nanmin(hiend_idx, axis=0)
                HIendCD = hiend_idx - pd + 1
                self.var._HIendCD[cond1] = HIendCD[cond1]

                # Duration of yield formation in calendar days
                self.var._YldFormCD[cond1] = (self.var._HIendCD - self.var._HIstartCD)[cond1]

                cond11 = (cond1 & (self.var._CropType == 3))

                # 1 Calendar days from sowing to end of flowering
                floweringend_idx = np.copy(day_idx)
                floweringend_idx[np.logical_not(GDDcum > self.var._FloweringEnd)] = 999
                # floweringend_idx[np.logical_not(GDDcum > self.var._FloweringEnd)] = np.nan
                floweringend_idx = np.nanmin(floweringend_idx, axis=0)
                FloweringEnd = floweringend_idx - pd + 1

                # 2 Duration of flowering in calendar days
                self.var._FloweringCD[cond11] = (FloweringEnd - self.var._HIstartCD)[cond11]

                # Harvest index growth coefficient
                self.calculate_HIGC()

                # Days to linear HI switch point
                self.calculate_HI_linear()

    def dynamic(self):
        """Function to update parameters for current crop grown as well 
        as counters pertaining to crop growth
        """
        # Update crop parameters for currently grown crops
        self.compute_water_productivity_adjustment_factor()
        self.adjust_planting_and_harvesting_date()
        self.update_growing_season()
        self.update_crop_parameters()
        self.select_crop_parameters()
        
class FAO56CropParameters(CropParameters):

    def initial(self):
        
        # Declare variables
        self.var.crop_parameters_to_read = [
            'PlantingDate','HarvestDate',
            'L_ini','L_dev','L_mid','L_late','Kc_ini','Kc_mid','Kc_end','p_std',#'Z',
            'Zmin','Zmax','Ky']

        self.var.crop_parameters_to_compute = [
            'L_ini_day','L_dev_day','L_mid_day','L_late_day',
            'PlantingDateAdj','HarvestDateAdj']

        # initialise parameters
        self.var.crop_parameter_names = self.var.crop_parameters_to_read + self.var.crop_parameters_to_compute
        arr_zeros = np.zeros((self.var.nCrop, self.var.nLat, self.var.nLon))
        for param in self.var.crop_parameter_names:
            nm = '_' + param
            vars(self.var)[nm] = arr_zeros

        # potential yield
        self.var.crop_parameter_names += ['Yx']
        self.var.PotYieldFileNC = self.var._configuration.cropOptions['PotYieldNC']
        self.var.PotYieldVarName = 'Yx' 
        if 'PotYieldVariableName' in self.var._configuration.cropOptions:
            self.var.PotYieldVarName = self.var._configuration.cropOptions['PotYieldVariableName']
        # self.var.co2_set_per_year  = False
        # self.var.Y = np.zeros((self.var.nRotation, self.var.nLat, self.var.nLon))
        
        # self.var.CropSequence = np.zeros((self.var.nCrop, self.var.nRotation, self.var.nLat, self.var.nLon))
        self.read()
        self.read_crop_sequence()

    def compute_growth_stage_length(self):
        nday = self.var._HarvestDateAdj - self.var._PlantingDateAdj
        self.var._L_ini_day = np.round(self.var._L_ini * nday)
        self.var._L_dev_day = np.round(self.var._L_dev * nday)
        self.var._L_mid_day = np.round(self.var._L_mid * nday)
        self.var._L_late_day = np.round(self.var._L_late * nday)  # TODO

    def read_potential_crop_yield(self):
        date = '%04i-%02i-%02i' %(self.var._modelTime.year, 1, 1)
        self.var._Yx = vos.netcdf2PCRobjClone(self.var.PotYieldFileNC,
                                              self.var.PotYieldVarName,
                                              date,
                                              useDoy = None,
                                              cloneMapFileName = self.var.cloneMap,
                                              LatitudeLongitude = True)
        
    def dynamic(self):
        self.adjust_planting_and_harvesting_date()
        self.compute_growth_stage_length()
        self.update_growing_season()
        self.read_potential_crop_yield()
        self.select_crop_parameters()
