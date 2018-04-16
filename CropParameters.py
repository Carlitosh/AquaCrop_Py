#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# AquaCrop
#
import os
import numpy as np
import pcraster as pcr
import virtualOS as vos
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

        self.cropParameterFileNC = iniItems.cropOptions['cropParameterNC']
        self.CalendarType = int(iniItems.cropOptions['CalendarType'])
        self.SwitchGDD = bool(int(iniItems.cropOptions['SwitchGDD']))
        self.GDDmethod = int(iniItems.cropOptions['GDDmethod'])

    def read_crop_parameters(self):
        """Function to read crop input parameters"""
        
        self.parameter_names = ['CropType','PlantingDate','HarvestDate','Emergence','MaxRooting','Senescence','Maturity','HIstart','Flowering','YldForm','PolHeatStress','PolColdStress','BioTempStress','PlantPop','Determinant','ETadj','LagAer','Tbase','Tupp','Tmax_up','Tmax_lo','Tmin_up','Tmin_lo','GDD_up','GDD_lo','fshape_b','PctZmin','Zmin','Zmax','fshape_r','fshape_ex','SxTopQ','SxBotQ','a_Tr','SeedSize','CCmin','CCx','CDC','CGC','Kcb','fage','WP','WPy','fsink','bsted','bface','HI0','HIini','dHI_pre','a_HI','b_HI','dHI0','exc','MaxFlowPct','p_up1','p_up2','p_up3','p_up4','p_lo1','p_lo2','p_lo3','p_lo4','fshape_w1','fshape_w2','fshape_w3','fshape_w4','Aer','beta','GermThr']

        for var in self.parameter_names:
            vars(self)[var] = vos.netcdf2PCRobjCloneWithoutTime(
                self.cropParameterFileNC,
                var,
                cloneMapFileName=self.cloneMap)

        self.nCrop = self.CropType.shape[0]

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

        # Update crop parameter names
        new_parameter_names = ['CC0','SxTop','SxBot']
        for nm in new_parameter_names:
            if nm not in self.parameter_names:
                self.parameter_names.append(nm)

    def compute_water_productivity_adjustment_factor(self, CO2):
        """Function to calculate water productivity adjustment factor 
        for elevation in C02 concentration"""

        # convenient to add crop dimension to CO2 variables
        CO2ref = CO2.ref[None,:,:] * np.ones((self.nCrop))[:,None,None]
        CO2conc = CO2.conc[None,:,:] * np.ones((self.nCrop))[:,None,None]
        
        # Get C02 weighting factor
        fw = np.zeros((self.nCrop, self.nLat, self.nLon))
        cond1 = (CO2conc > CO2ref)
        cond11 = (cond1 & (CO2conc >= 550))
        fw[cond11] = 1
        cond12 = (cond1 & np.logical_not(cond11))
        fw[cond12] = 1 - ((550 - CO2conc) / (550 - CO2ref))

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
        
    # def growing_degree_day(self, currTimeStep, meteo):
    #     """Function to calculate number of growing degree days on 
    #     current day
    #     """
    #     tmax = meteo.tmax[None,:,:] * np.ones((self.nCrop))[:,None,None]
    #     tmin = meteo.tmin[None,:,:] * np.ones((self.nCrop))[:,None,None]
        
    #     if self.GDDmethod == 1:
    #         tmean = ((tmax + tmin) / 2)
    #         tmean = np.clip(tmean, self.Tbase, self.Tupp)
    #     elif self.GDDmethod == 2:
    #         tmax = np.clip(tmax, self.Tbase, self.Tupp)
    #         tmin = np.clip(tmin, self.Tbase, self.Tupp)
    #         tmean = ((tmax + tmin) / 2)
    #     elif self.GDDmethod == 3:
    #         tmax = np.clip(tmax, self.Tbase, self.Tupp)
    #         tmin = np.clip(tmin, None, self.Tupp)
    #         tmean = np.clip(tmean, self.Tbase, None)
            
    #     self.GDD = (tmean - self.Tbase)
    #     # TODO: why have you commented this out?
    #     # self.GDDcum += GDD

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

    def update_growing_degree_day(self, Meteo, startTime, MASK=None):
        """Function to compute growing degree day parameter for upcoming 
        growing season"""

        if MASK is None:
            MASK = np.full(self.PlantingDate.shape, True)

        if not np.any(MASK):
            GDD = np.full(self.PlantingDate.shape, np.nan)
        else:
            st = startTime
            sd = startTime.timetuple().tm_yday # (Julian day: 1 Jan = 1)

            # objective of the following code section is to obtain an
            # index of the first complete growing season of the given
            # crop in each grid cell

            # planting/harvesting date of crop
            pd = self.PlantingDate
            hd = self.HarvestDate

            # if start day of simulation is greater than planting day the
            # first complete growing season will not be until the
            # following year
            hd[sd > pd] += 365
            pd[sd > pd] += 365

            # if start day is less than or equal to planting day, but
            # harvest day is less than planting day, the harvest day will
            # occur in the following year
            hd[((sd <= pd) & (hd < pd))] += 365

            # adjust values for leap year
            isLeapYear1 = calendar.isleap(st.year)
            isLeapYear2 = calendar.isleap(st.year + 1)
            pd[(isLeapYear1 & (pd >= 60))] += 1  # TODO: check these
            hd[(isLeapYear1 & (hd >= 60))] += 1
            pd[(isLeapYear2 & (pd >= 425))] += 1
            hd[(isLeapYear2 & (hd >= 425))] += 1

            pd[np.logical_not(MASK)] = np.nan
            hd[np.logical_not(MASK)] = np.nan

            ndays = np.nanmax(hd)  # numpy.nanmax calculates the maximum value, ignoring NaN cell values
            day_idx = np.arange(1, ndays + 1)[:,None,None,None] * np.ones((self.nCrop, self.nLon, self.nLat))[None,:,:,:]

            # Dimensions of pd are crop,lat,lon; day_idx dims are
            # time,crop,lat,lon - broadcasting should work automatically
            growing_season_idx = ((pd >= day_idx) & (hd <= day_idx))

            # Extract weather data for first growing season
            tmin = vos.netcdf2NumPyTimeSlice(meteo.tmpFileNC, meteo.tmnVarName,
                                             startTime,
                                             startTime + datetime.timedelta(ndays),
                                             cloneMapFileName = self.cloneMap,
                                             LatitudeLongitude = True)

            tmax = vos.netcdf2NumPyTimeSlice(meteo.tmpFileNC, meteo.tmxVarName,
                                             startTime,
                                             startTime + datetime.timedelta(ndays),
                                             cloneMapFileName = self.cloneMap,
                                             LatitudeLongitude = True)

            # broadcast to crop dimension
            tmax = tmax[:,None,:,:] * np.ones((self.nCrop))[None,:,None,None]
            tmin = tmin[:,None,:,:] * np.ones((self.nCrop))[None,:,None,None] 

            # for convenience
            tupp = self.Tupp[None,:,:,:] * np.ones((ndays))[:,None,None,None]
            tbase = self.Tbase[None,:,:,:] * np.ones((ndays))[:,None,None,None]

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
                tmean = np.clip(tmean, self.Tbase, None)

            tmean[np.logical_not(growing_season_idx)] = 0
            tbase[np.logical_not(growing_season_idx)] = 0
            GDD = (tmean - tbase)

        return GDD

    
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
        self.FloweringEnd = np.zeros((self.nCrop, self.nLat, self.nLon))  # TODO: check not used unless CropType == 3
        self.FloweringEnd[cond2] = (self.HIstart + self.Flowering)[cond2]
        
        Mode = self.CalendarType
        if Mode == 1:

            # "Duplicate calendar values (needed to minimise if-statements when
            # switching between GDD and CD runs)
            self.EmergenceCD = self.Emergence
            self.Canopy10PctCD = self.Canopy10Pct
            self.MaxRootingCD = self.MaxRooting
            self.SenescenceCD = self.Senescence
            self.MaturityCD = self.Maturity
            self.MaxCanopyCD = self.MaxCanopy
            self.CanopyDevEndCD = self.CanopyDevEnd
            self.HIstartCD = self.HIstart
            self.HIendCD = self.Hiend
            self.YldFormCD = self.YldForm            
            self.FloweringEndCD = self.FloweringEnd
            self.FloweringCD = self.Flowering

        # Pre-compute cumulative GDD during growing season
        if (Mode == 1 & self.SwitchGDD) | (Mode == 2):
            GDD = update_growing_degree_day(self, Meteo, currTimeStep.startTime)
            GDDcum = np.cumsum(GDD, axis=0)

            # "Check if converting crop calendar to GDD mode"
            if Mode == 1 & self.SwitchGDD:

                # Find GDD equivalent for each crop calendar variable
                m,n,p = pd.shape  # crop,lat,lon
                I,J,K = np.ogrid[:m,:n,:p]

                emergence_idx = pd + self.EmergenceCD  # crop,lat,lon
                self.Emergence = GDDcum[emergence_idx,I,J,K]
                canopy10pct_idx = pd + self.Canopy10pctCD
                self.Canopy10pct = GDDcum[canopy10pct_idx,I,J,K]
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
                self.CGC_CD = self.CGC
                self.CGC = (np.log((((0.98 * self.CCx) - self.CCx) * self.CC0) / (-0.25 * (self.CCx ** 2)))) / (-(self.MaxCanopy - self.Emergence))

                # "Convert CDC to GDD mode"
                self.CDC_CD = self.CDC
                tCD = self.MaturityCD - self.SenescenceCD
                tCD[tCD <= 0] = 1
                tGDD = self.Maturity - self.Senescence
                tGDD[tGDD <= 0] = 5
                self.CDC = (self.CCx / tGDD) * np.log(1 + ((1 - self.CCi / self.CCx) / 0.05))

                # "Set calendar type to GDD mode"
                iniItems.cropOptions['CalendarType'] = "2"
            
            elif Mode == 2:

                # "Find calendar days [equivalent] for some variables"

                # "1 Calendar days from sowing to maximum canopy cover"
                maxcanopy_idx = day_idx
                maxcanopy_idx[np.logical_not(GDDcum > self.MaxCanopy)] = np.nan
                maxcanopy_idx = np.nanmin(maxcanopy_idx, axis=0)
                self.MaxCanopyCD = maxcanopy_idx - pd + 1

                # "2 Calendar days from sowing to end of vegetative growth"
                canopydevend_idx = day_idx
                canopydevend_idx[np.logical_not(GDDcum > self.CanopyDevEnd)] = np.nan
                canopydevend_idx = np.nanmin(canopydevend_idx, axis=0)
                self.CanopyDevEndCD = canopydevend_idx - pd + 1

                # "3 Calendar days from sowing to start of yield formation"
                histart_idx = day_idx
                histart_idx[np.logical_not(GDDcum > self.HIstart)] = np.nan
                histart_idx = np.nanmin(histart_idx, axis=0)
                self.HIstartCD = histart_idx - pd + 1

                # "4 Calendar days from sowing to end of yield formation"
                hiend_idx = day_idx
                hiend_idx[np.logical_not(GDDcum > self.HIend)] = np.nan
                hiend_idx = np.nanmin(hiend_idx, axis=0)
                self.HIendCD = hiend_idx - pd + 1

                # "Duration of yield formation in calendar days"
                self.YldFormCD = self.HIendCD - self.HIstartCD

                cond1 = (self.CropType == 3)

                # "1 Calendar days from sowing to end of flowering"
                floweringend_idx = day_idx
                floweringend_idx[np.logical_not(GDDcum > self.FloweringEnd)] = np.nan
                floweringend_idx = np.nanmin(floweringend_idx, axis=0)
                FloweringEnd = floweringend_idx - pd + 1

                # "2 Duration of flowering in calendar days"
                self.FloweringCD[cond1] = (FloweringEnd - self.HIstartCD)[cond1]

    def update_crop_parameters(self, currTimeStep, Meteo):

        # Identify crops which are growing in the cell and planted on the
        # current day
        cond1 = (self.CropSequence & (currTimeStep.doy == self.PlantingDate))

        GDD = update_growing_degree_day(self, Meteo, currTimeStep.fulldate, MASK=cond1)
        GDDcum = np.cumsum(GDD, axis=0)

        day_idx = np.arange(1, GDD.shape[0])[:,None,None,None] * np.ones((self.nCrop, self.nLon, self.nLat))[None,:,:,:]
        
        # Find calendar days [equivalent] for some variables

        # 1 - Calendar days from sowing to maximum canopy cover
        maxcanopy_idx = day_idx
        maxcanopy_idx[np.logical_not(GDDcum > self.MaxCanopy)] = np.nan
        maxcanopy_idx = np.nanmin(maxcanopy_idx, axis=0)
        MaxCanopyCD = (maxcanopy_idx - pd + 1)
        self.MaxCanopyCD[cond1] = MaxCanopyCD[cond1]

        # 2 - Calendar days from sowing to end of vegetative growth
        canopydevend_idx = day_idx
        canopydevend_idx[np.logical_not(GDDcum > self.CanopyDevEnd)] = np.nan
        canopydevend_idx = np.nanmin(canopydevend_idx, axis=0)
        CanopyDevEndCD = canopydevend_idx - pd + 1
        self.CanopyDevEndCD[cond1] = CanopyDevEndCD[cond1]

        # 3 - Calendar days from sowing to start of yield formation
        histart_idx = day_idx
        histart_idx[np.logical_not(GDDcum > self.HIstart)] = np.nan
        histart_idx = np.nanmin(histart_idx, axis=0)
        HIstartCD = histart_idx - pd + 1
        self.HIstartCD[cond1] = HIstartCD[cond1]
        
        # 4 - Calendar days from sowing to end of yield formation
        hiend_idx = day_idx
        hiend_idx[np.logical_not(GDDcum > self.HIend)] = np.nan
        hiend_idx = np.nanmin(hiend_idx, axis=0)
        HIendCD = hiend_idx - pd + 1
        self.HIendCD[cond1] = HIendCD[cond1]
        
        # Duration of yield formation in calendar days
        self.YldFormCD[cond1] = (self.HIendCD - self.HIstartCD)[cond1]

        cond11 = (cond1 & (self.CropType == 3))

        # 1 Calendar days from sowing to end of flowering
        floweringend_idx = day_idx
        floweringend_idx[np.logical_not(GDDcum > self.FloweringEnd)] = np.nan
        floweringend_idx = np.nanmin(floweringend_idx, axis=0)
        FloweringEnd = floweringend_idx - pd + 1

        # 2 Duration of flowering in calendar days
        self.FloweringCD[cond11] = (FloweringEnd - self.HIstartCD)[cond11]

        # Harvest index growth coefficient
        self.AOS_CalculateHIGC()

        # Days to linear HI switch point
        self.AOS_CalculateHILinear()
                
