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

    def __init__(self, iniItems, landmask):
        object.__init__(self)

        self.cloneMap = iniItems.cloneMap
        self.landmask = landmask

        attr = vos.getMapAttributesALL(self.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])

        self.nCrop = int(iniItems.cropOptions['nCrop'])
        self.nRotation = int(iniItems.cropOptions['nRotation'])
        self.cropParameterFileNC = str(iniItems.cropOptions['cropParameterNC'])
        self.CalendarType = int(iniItems.cropOptions['CalendarType'])
        self.SwitchGDD = bool(int(iniItems.cropOptions['SwitchGDD']))
        self.GDDmethod = int(iniItems.cropOptions['GDDmethod'])
        
        # Declare variables
        self.crop_parameters_to_read = [
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
        
        self.crop_parameters_to_compute = [
            'CC0','SxTop','SxBot','fCO2','tLinSwitch','dHILinear','HIGC',
            'CanopyDevEnd','CanopyDevEndCD','Canopy10Pct','MaxCanopy',
            'MaxCanopyCD','HIstartCD','HIend','HIendCD','YldFormCD',
            'FloweringEnd','Flowering','FloweringCD','CGC','CDC']

        self.crop_parameter_names = self.crop_parameters_to_read + self.crop_parameters_to_compute
        arr_zeros = np.zeros((self.nCrop, self.nLat, self.nLon))
        for param in self.crop_parameter_names:
            vars(self)[param] = arr_zeros

        self.PlantingDateAdj = arr_zeros
        self.HarvestDateAdj = arr_zeros
        self.CropSequence = np.zeros((self.nCrop, self.nRotation, self.nLat, self.nLon))
        
    def read(self):
        """Function to read crop input parameters"""
        
        for param in self.crop_parameters_to_read:
            vars(self)[param] = vos.netcdf2PCRobjCloneWithoutTime(
                self.cropParameterFileNC,
                param,
                cloneMapFileName=self.cloneMap)

        # self.nCrop = self.CropType.shape[0]

        # Crop sequence - i.e. the order in which crops are grown
        CropSequence = vos.netcdf2PCRobjCloneWithoutTime(
            self.cropParameterFileNC,
            'CropSequence',
            cloneMapFileName=self.cloneMap)

        CropSequence = (
            CropSequence[:,:,None,None]
            * np.ones((self.nLat, self.nLon))[None,None,:,:])
        self.CropSequence = CropSequence.astype(bool)

        # self.nCrop = self.CropSequence.shape[0]
        # self.nRotation = self.CropSequence.shape[1]
        
    def compute_variables(self, currTimeStep, meteo):
        """Function to compute additional crop variables required to run 
        AquaCrop"""
        
        # The following adapted from lines 160-240 of AOS_ComputeVariables.m

        # Fractional canopy cover size at emergence
        self.CC0 = np.round(10000 * self.PlantPop * self.SeedSize * 10 ** -8) / 10000
        
        # Root extraction terms
        S1 = self.SxTopQ
        S2 = self.SxBotQ
        self.SxTop = np.zeros((self.nCrop, self.nLat, self.nLon))
        self.SxBot = np.zeros((self.nCrop, self.nLat, self.nLon))

        cond1 = (S1 == S2)
        self.SxTop[cond1] = S1[cond1]
        self.SxBot[cond1] = S2[cond1]

        cond2 = np.logical_not(cond1)
        cond21 = (cond2 & (self.SxTopQ < self.SxBotQ))
        S1[cond21] = self.SxBotQ[cond21]
        S2[cond21] = self.SxTopQ[cond21]
        xx = 3 * (S2 / (S1 - S2))
        SS1 = np.zeros((self.nCrop, self.nLat, self.nLon))
        SS2 = np.zeros((self.nCrop, self.nLat, self.nLon))

        cond22 = (cond2 & (xx < 0.5))
        SS1[cond22] = ((4 / 3.5) * S1)[cond22]
        SS2[cond22] = 0

        cond23 = (cond2 & np.logical_not(cond22))
        SS1[cond23] = ((xx + 3.5) * (S1 / (xx + 3)))[cond23]
        SS2[cond23] = ((xx - 0.5) * (S2 / xx))[cond23]

        cond24 = (cond2 & (self.SxTopQ > self.SxBotQ))
        self.SxTop[cond24] = SS1[cond24]
        self.SxBot[cond24] = SS2[cond24]

        cond25 = (cond2 & np.logical_not(cond24))
        self.SxTop[cond25] = SS2[cond25]
        self.SxBot[cond25] = SS1[cond25]

        # Crop calender
        self.AOS_ComputeCropCalendar(currTimeStep, meteo)

        # Harvest index growth coefficient
        self.AOS_CalculateHIGC()

        # Days to linear HI switch point
        self.AOS_CalculateHILinear()

        # # Update crop parameter names
        # new_parameter_names = ['CC0','SxTop','SxBot']
        # for nm in new_parameter_names:
        #     if nm not in self.parameter_names:
        #         self.parameter_names.append(nm)

    def compute_water_productivity_adjustment_factor(self, CO2):
        """Function to calculate water productivity adjustment factor 
        for elevation in CO2 concentration"""

        # convenient to add crop dimension to CO2 variables
        CO2ref = CO2.RefConc
        CO2conc = CO2.conc[None,:,:] * np.ones((self.nCrop))[:,None,None]
        
        # Get CO2 weighting factor
        fw = np.zeros((self.nCrop, self.nLat, self.nLon))
        cond1 = (CO2conc > CO2ref)
        cond11 = (cond1 & (CO2conc >= 550))
        fw[cond11] = 1
        cond12 = (cond1 & np.logical_not(cond11))
        fw[cond12] = (1 - ((550 - CO2conc) / (550 - CO2ref)))[cond12]

        # Determine adjustment for each crop in first year of simulation
        fCO2 = ((CO2conc / CO2ref) /
                (1 + (CO2conc - CO2ref) * ((1 - fw)
                                           * self.bsted + fw
                                           * ((self.bsted * self.fsink)
                                              + (self.bface
                                                 * (1 - self.fsink))))))

        # Consider crop type
        ftype = (40 - self.WP) / (40 - 20)
        ftype = np.clip(ftype, 0, 1)
        self.fCO2 = 1 + ftype * (fCO2 - 1)
        
    def AOS_CalculateHILinear(self):
        """Function to calculate time to switch to linear harvest index 
        build-up, and associated linear rate of build-up. Only for 
        fruit/grain crops
        """
        # Determine linear switch point
        cond1 = (self.CropType == 3)
        ti = np.zeros((self.nCrop, self.nLat, self.nLon))
        tmax = self.YldFormCD
        HIest = np.zeros((self.nCrop, self.nLat, self.nLon))
        HIprev = self.HIini

        # Iterate to find linear switch point
        while np.any((self.CropType == 3) & (HIest <= self.HI0) & (ti < tmax)):    
            ti += 1
            HInew = ((self.HIini * self.HI0) / (self.HIini + (self.HI0 - self.HIini) * np.exp(-self.HIGC * ti)))
            HIest = (HInew + (tmax - ti) * (HInew - HIprev))
            HIprev = HInew

        self.tLinSwitch = ti
        self.tLinSwitch[self.CropType != 3] = np.nan
            
        # Determine linear build-up rate
        cond1 = (self.tLinSwitch > 0)
        HIest[cond1] = ((self.HIini * self.HI0) / (self.HIini + (self.HI0 - self.HIini) * np.exp(-self.HIGC * self.tLinSwitch)))[cond1]
        HIest[np.logical_not(cond1)] = 0
        self.dHILinear = ((self.HI0 - HIest) / (tmax - self.tLinSwitch))  # dHILin will be set to nan in the same cells as tSwitch
                    
    def AOS_CalculateHIGC(self):
        """Function to calculate harvest index growth coefficient"""
        # Total yield formation days
        tHI = self.YldFormCD

        # Iteratively estimate HIGC
        self.HIGC = np.full((self.nCrop, self.nLat, self.nLon), 0.001)
        HIest = np.zeros((self.nCrop, self.nLat, self.nLon))
        while np.any(HIest < (0.98 * self.HI0)):
            cond1 = (HIest < (0.98 * self.HI0))
            self.HIGC += 0.001
            HIest = ((self.HIini * self.HI0) / (self.HIini + (self.HI0 - self.HIini) * np.exp(-self.HIGC * tHI)))
            
        self.HIGC[HIest >= self.HI0] -= 0.001

    def AOS_ComputeCropCalendar(self, currTimeStep, meteo):
       
        # "Time from sowing to end of vegetative growth period"
        cond1 = (self.Determinant == 1)
        self.CanopyDevEnd = self.Senescence
        self.CanopyDevEnd[cond1] = (np.round(self.HIstart + (self.Flowering / 2)))[cond1]

        # "Time from sowing to 10% canopy cover (non-stressed conditions)
        self.Canopy10Pct = np.round(self.Emergence + (np.log(0.1 / self.CC0) / self.CGC))

        # "Time from sowing to maximum canopy cover (non-stressed conditions)
        self.MaxCanopy = np.round(self.Emergence + (np.log((0.25 * self.CCx * self.CCx / self.CC0) / (self.CCx - (0.98 * self.CCx))) / self.CGC))

        # "Time from sowing to end of yield formation"
        self.HIend = self.HIstart + self.YldForm

        cond2 = (self.CropType == 3)

        # TODO: declare these in __init__
        self.FloweringEnd = np.zeros((self.nCrop, self.nLat, self.nLon))
        self.FloweringEndCD = np.zeros((self.nCrop, self.nLat, self.nLon))
        self.FloweringCD = np.zeros((self.nCrop, self.nLat, self.nLon))
                
        self.FloweringEnd[cond2] = (self.HIstart + self.Flowering)[cond2]
        self.FloweringEndCD[cond2] = self.FloweringEndCD[cond2]
        self.FloweringCD[cond2] = self.Flowering[cond2]
        
        # Mode = self.CalendarType
        # if Mode == 1:
        if self.CalendarType == 1:
            # "Duplicate calendar values (needed to minimise if-statements when
            # switching between GDD and CD runs)
            self.EmergenceCD = np.copy(self.Emergence)      # only used in this function
            self.Canopy10PctCD = np.copy(self.Canopy10Pct)  # only used in this function
            self.MaxRootingCD = np.copy(self.MaxRooting)    # only used in this function
            self.SenescenceCD = np.copy(self.Senescence)    # only used in this function
            self.MaturityCD = np.copy(self.Maturity)        # only used in this function
            self.MaxCanopyCD = np.copy(self.MaxCanopy)
            self.CanopyDevEndCD = np.copy(self.CanopyDevEnd)
            self.HIstartCD = np.copy(self.HIstart)
            self.HIendCD = np.copy(self.HIend)
            self.YldFormCD = np.copy(self.YldForm)            
            self.FloweringEndCD = np.copy(self.FloweringEnd)  # only used in this function
            self.FloweringCD = np.copy(self.Flowering)

        # Pre-compute cumulative GDD during growing season
        if (self.CalendarType == 1 & self.SwitchGDD) | (self.CalendarType == 2):

            pd = np.copy(self.PlantingDate)
            hd = np.copy(self.HarvestDate)
            sd = currTimeStep.startTime.timetuple().tm_yday

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
            isLeapYear1 = calendar.isleap(currTimeStep.startTime.year)
            isLeapYear2 = calendar.isleap(currTimeStep.startTime.year + 1)
            pd[(isLeapYear1 & (pd >= 60))] += 1  # TODO: check these
            hd[(isLeapYear1 & (hd >= 60))] += 1
            pd[(isLeapYear2 & (pd >= 425))] += 1
            hd[(isLeapYear2 & (hd >= 425))] += 1            

            max_harvest_date = int(np.max(hd))
            day_idx = np.arange(sd, max_harvest_date + 1)[:,None,None,None] * np.ones((self.nCrop, self.nLon, self.nLat))[None,:,:,:]
            growing_season_idx = ((day_idx >= pd) & (day_idx <= hd))

            # Extract weather data for first growing season
            tmin = vos.netcdf2NumPyTimeSlice(meteo.tmpFileNC, meteo.tmnVarName,
                                             currTimeStep.startTime,
                                             currTimeStep.startTime + datetime.timedelta(int(max_harvest_date - sd)),
                                             cloneMapFileName = self.cloneMap,
                                             LatitudeLongitude = True)

            tmax = vos.netcdf2NumPyTimeSlice(meteo.tmpFileNC, meteo.tmxVarName,
                                             currTimeStep.startTime,
                                             currTimeStep.startTime + datetime.timedelta(int(max_harvest_date - sd)),
                                             cloneMapFileName = self.cloneMap,
                                             LatitudeLongitude = True)

            # broadcast to crop dimension
            tmax = tmax[:,None,:,:] * np.ones((self.nCrop))[None,:,None,None]
            tmin = tmin[:,None,:,:] * np.ones((self.nCrop))[None,:,None,None]

            # for convenience
            tupp = self.Tupp[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]
            tbase = self.Tbase[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]

            # calculate GDD according to the various methods
            if self.GDDmethod == 1:
                tmean = ((tmax + tmin) / 2)
                tmean = np.clip(tmean, self.Tbase, self.Tupp)
            elif self.GDDmethod == 2:
                tmax = np.clip(tmax, self.Tbase, self.Tupp)
                tmin = np.clip(tmin, self.Tbase, self.Tupp)
                tmean = ((tmax + tmin) / 2)
            elif self.GDDmethod == 3:
                tmax = np.clip(tmax, self.Tbase, self.Tupp)
                tmin = np.clip(tmin, None, self.Tupp)
                tmean = ((tmax + tmin) / 2)
                tmean = np.clip(tmean, self.Tbase, None)

            tmean[np.logical_not(growing_season_idx)] = 0
            tbase[np.logical_not(growing_season_idx)] = 0
            GDD = (tmean - tbase)
            GDDcum = np.cumsum(GDD, axis=0)

            # "Check if converting crop calendar to GDD mode"
            # if Mode == 1 & self.SwitchGDD:
            if self.CalendarType == 1 & self.SwitchGDD:

                # Find GDD equivalent for each crop calendar variable
                m,n,p = pd.shape  # crop,lat,lon
                I,J,K = np.ogrid[:m,:n,:p]

                emergence_idx = pd + self.EmergenceCD  # crop,lat,lon
                self.Emergence = GDDcum[emergence_idx,I,J,K]
                canopy10pct_idx = pd + self.Canopy10PctCD
                self.Canopy10Pct = GDDcum[canopy10pct_idx,I,J,K]
                maxrooting_idx = pd + self.MaxRootingCD
                self.MaxRooting = GDDcum[maxrooting_idx,I,J,K]
                maxcanopy_idx = pd + self.MaxCanopyCD
                self.MaxCanopy = GDDcum[maxcanopy_idx,I,J,K]
                canopydevend_idx = pd + self.CanopyDevEndCD
                self.CanopyDevEnd = GDDcum[canopydevend_idx,I,J,K]
                senescence_idx = pd + self.SenescenceCD
                self.Senescence = GDDcum[senescence_idx,I,J,K]
                maturity_idx = pd + self.MaturityCD
                self.Maturity = GDDcum[maturity_idx,I,J,K]
                histart_idx = pd + self.HIstartCD
                self.HIstart = GDDcum[histart_idx,I,J,K]
                hiend_idx = pd + self.HIendCD
                self.HIend = GDDcum[hiend_idx,I,J,K]
                yldform_idx = pd + self.YldFormCD
                self.YldForm = GDDcum[yldform_idx,I,J,K]

                cond2 = (self.CropType == 3)
                floweringend_idx = pd + self.FloweringEndCD
                self.FloweringEnd[cond2] = GDDcum[floweringend_idx,I,J,K][cond2]
                self.Flowering[cond2] = (self.FloweringEnd - self.HIstart)[cond2]

                # "Convert CGC to GDD mode"
                # self.CGC_CD = self.CGC
                self.CGC = (np.log((((0.98 * self.CCx) - self.CCx) * self.CC0) / (-0.25 * (self.CCx ** 2)))) / (-(self.MaxCanopy - self.Emergence))

                # "Convert CDC to GDD mode"
                # self.CDC_CD = self.CDC
                tCD = self.MaturityCD - self.SenescenceCD
                tCD[tCD <= 0] = 1
                tGDD = self.Maturity - self.Senescence
                tGDD[tGDD <= 0] = 5
                self.CDC = (self.CCx / tGDD) * np.log(1 + ((1 - self.CCi / self.CCx) / 0.05))

                # "Set calendar type to GDD mode"
                iniItems.cropOptions['CalendarType'] = "2"
            
            # elif Mode == 2:
            elif self.CalendarType == 2:

                # "Find calendar days [equivalent] for some variables"

                # "1 Calendar days from sowing to maximum canopy cover"

                # # TODO: check this indexing
                # day_idx = np.arange(0, GDD.shape[0])[:,None,None,None] * np.ones((self.nCrop, self.nLon, self.nLat))[None,:,:,:]
                # pd,hd = self.adjust_planting_and_harvesting_date(currTimeStep.startTime)

                maxcanopy_idx = np.copy(day_idx)
                maxcanopy_idx[np.logical_not(GDDcum > self.MaxCanopy)] = np.nan
                maxcanopy_idx = np.nanmin(maxcanopy_idx, axis=0)
                self.MaxCanopyCD = maxcanopy_idx - pd + 1

                # "2 Calendar days from sowing to end of vegetative growth"
                canopydevend_idx = np.copy(day_idx)
                canopydevend_idx[np.logical_not(GDDcum > self.CanopyDevEnd)] = np.nan
                canopydevend_idx = np.nanmin(canopydevend_idx, axis=0)
                self.CanopyDevEndCD = canopydevend_idx - pd + 1

                # "3 Calendar days from sowing to start of yield formation"
                histart_idx = np.copy(day_idx)
                histart_idx[np.logical_not(GDDcum > self.HIstart)] = np.nan
                histart_idx = np.nanmin(histart_idx, axis=0)
                self.HIstartCD = histart_idx - pd + 1

                # "4 Calendar days from sowing to end of yield formation"
                hiend_idx = np.copy(day_idx)
                hiend_idx[np.logical_not(GDDcum > self.HIend)] = np.nan
                hiend_idx = np.nanmin(hiend_idx, axis=0)
                self.HIendCD = hiend_idx - pd + 1

                # "Duration of yield formation in calendar days"
                self.YldFormCD = self.HIendCD - self.HIstartCD

                cond1 = (self.CropType == 3)

                # "1 Calendar days from sowing to end of flowering"
                floweringend_idx = np.copy(day_idx)
                floweringend_idx[np.logical_not(GDDcum > self.FloweringEnd)] = np.nan
                floweringend_idx = np.nanmin(floweringend_idx, axis=0)
                FloweringEnd = floweringend_idx - pd + 1

                # "2 Duration of flowering in calendar days"
                self.FloweringCD[cond1] = (FloweringEnd - self.HIstartCD)[cond1]

    def update(self, currTimeStep, meteo):
        """Function to update certain crop parameters for current 
        time step
        """
        pd = np.copy(self.PlantingDate)
        hd = np.copy(self.HarvestDate)
        st = currTimeStep.currTime
        sd = currTimeStep.currTime.timetuple().tm_yday

        # adjust values for leap year (objective is to preserve date)
        isLeapYear1 = calendar.isleap(st.year)
        isLeapYear2 = calendar.isleap(st.year + 1)
        pd[(isLeapYear1 & (pd >= 60))] += 1  # TODO: check these
        hd[(isLeapYear1 & (hd >= 60) & (hd > pd))] += 1
        hd[(isLeapYear2 & (hd >= 60) & (hd < pd))] += 1

        cond1 = ((pd <= currTimeStep.doy) & (currTimeStep.doy <= hd))
        cond2 = ((pd > hd) & ((pd <= currTimeStep.doy) | (currTimeStep.doy <= hd)))

        # if harvest day is less than planting day, harvest will take place
        # in following year (leap year already accounted for, so add 365)
        hd[hd < pd] += 365
        hd_arr = (np.datetime64(str(datetime.datetime(currTimeStep.year, 1, 1))) + np.array(hd - 1, dtype='timedelta64[D]'))
        cond3 = hd_arr <= np.datetime64(str(currTimeStep.endTime))        
        # td = np.array(hd - pd, dtype='timedelta64[D]')
        # cond3 = ((np.datetime64(str(currTimeStep.currTime)) + td) <= np.datetime64(str(currTimeStep.endTime)))
        self.GrowingSeason = ((cond1 & cond3) | (cond2 & cond3))
        
        # Update certain crop parameters if using GDD mode
        if (self.CalendarType == 2):

            cond1 = ((self.GrowingSeason) & (currTimeStep.doy == pd))
            pd[np.logical_not(cond1)] = 0
            hd[np.logical_not(cond1)] = 0
            max_harvest_date = int(np.max(hd))

            if (max_harvest_date > 0):
                
                day_idx = np.arange(sd, max_harvest_date + 1)[:,None,None,None] * np.ones((self.nCrop, self.nLon, self.nLat))[None,:,:,:]
                growing_season_idx = ((day_idx >= pd) & (day_idx <= hd))

                # Extract weather data for first growing season
                tmin = vos.netcdf2NumPyTimeSlice(meteo.tmpFileNC, meteo.tmnVarName,
                                                 currTimeStep.currTime,
                                                 currTimeStep.currTime + datetime.timedelta(int(max_harvest_date - sd)),
                                                 cloneMapFileName = self.cloneMap,
                                                 LatitudeLongitude = True)

                tmax = vos.netcdf2NumPyTimeSlice(meteo.tmpFileNC, meteo.tmxVarName,
                                                 currTimeStep.currTime,
                                                 currTimeStep.currTime + datetime.timedelta(int(max_harvest_date - sd)),
                                                 cloneMapFileName = self.cloneMap,
                                                 LatitudeLongitude = True)

                # broadcast to crop dimension
                tmax = tmax[:,None,:,:] * np.ones((self.nCrop))[None,:,None,None]
                tmin = tmin[:,None,:,:] * np.ones((self.nCrop))[None,:,None,None]

                # for convenience
                tupp = self.Tupp[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]
                tbase = self.Tbase[None,:,:,:] * np.ones((tmin.shape[0]))[:,None,None,None]

                # calculate GDD according to the various methods
                if self.GDDmethod == 1:
                    tmean = ((tmax + tmin) / 2)
                    tmean = np.clip(tmean, self.Tbase, self.Tupp)
                elif self.GDDmethod == 2:
                    tmax = np.clip(tmax, self.Tbase, self.Tupp)
                    tmin = np.clip(tmin, self.Tbase, self.Tupp)
                    tmean = ((tmax + tmin) / 2)
                elif self.GDDmethod == 3:
                    tmax = np.clip(tmax, self.Tbase, self.Tupp)
                    tmin = np.clip(tmin, None, self.Tupp)
                    tmean = ((tmax + tmin) / 2)
                    tmean = np.clip(tmean, self.Tbase, None)

                tmean[np.logical_not(growing_season_idx)] = 0
                tbase[np.logical_not(growing_season_idx)] = 0
                GDD = (tmean - tbase)
                GDDcum = np.cumsum(GDD, axis=0)
        
                # 1 - Calendar days from sowing to maximum canopy cover
                maxcanopy_idx = np.copy(day_idx)
                maxcanopy_idx[np.logical_not(GDDcum > self.MaxCanopy)] = np.nan
                maxcanopy_idx = np.nanmin(maxcanopy_idx, axis=0)
                MaxCanopyCD = (maxcanopy_idx - pd + 1)
                self.MaxCanopyCD[cond1] = MaxCanopyCD[cond1]

                # 2 - Calendar days from sowing to end of vegetative growth
                canopydevend_idx = np.copy(day_idx)
                canopydevend_idx[np.logical_not(GDDcum > self.CanopyDevEnd)] = np.nan
                canopydevend_idx = np.nanmin(canopydevend_idx, axis=0)
                CanopyDevEndCD = canopydevend_idx - pd + 1
                self.CanopyDevEndCD[cond1] = CanopyDevEndCD[cond1]

                # 3 - Calendar days from sowing to start of yield formation
                histart_idx = np.copy(day_idx)
                histart_idx[np.logical_not(GDDcum > self.HIstart)] = np.nan
                histart_idx = np.nanmin(histart_idx, axis=0)
                HIstartCD = histart_idx - pd + 1
                self.HIstartCD[cond1] = HIstartCD[cond1]

                # 4 - Calendar days from sowing to end of yield formation
                hiend_idx = np.copy(day_idx)
                hiend_idx[np.logical_not(GDDcum > self.HIend)] = np.nan
                hiend_idx = np.nanmin(hiend_idx, axis=0)
                HIendCD = hiend_idx - pd + 1
                self.HIendCD[cond1] = HIendCD[cond1]

                # Duration of yield formation in calendar days
                self.YldFormCD[cond1] = (self.HIendCD - self.HIstartCD)[cond1]

                cond11 = (cond1 & (self.CropType == 3))

                # 1 Calendar days from sowing to end of flowering
                floweringend_idx = np.copy(day_idx)
                floweringend_idx[np.logical_not(GDDcum > self.FloweringEnd)] = np.nan
                floweringend_idx = np.nanmin(floweringend_idx, axis=0)
                FloweringEnd = floweringend_idx - pd + 1

                # 2 Duration of flowering in calendar days
                self.FloweringCD[cond11] = (FloweringEnd - self.HIstartCD)[cond11]

                # Harvest index growth coefficient
                self.AOS_CalculateHIGC()

                # Days to linear HI switch point
                self.AOS_CalculateHILinear()
                
