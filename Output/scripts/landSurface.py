#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# AquaCrop
#

import types
import pcraster as pcr
import virtualOS as vos
import parameterSoilAndTopo as parSoilAndTopo
import parameterCrop as parCrop

import logging
logger = logging.getLogger(__name__)

from ncConverter import *

class LandSurface(object):
    
    def getState(self):
        # TODO: use PCRGLOBWB as example
        pass
    
    def getPseudoState(self):
        # TODO: use PCRGLOBWB as example
        pass

    def __init__(self,iniItems,currTimeStep,meteo,groundwater,landmask,initialState=None):

        object.__init__(self)

        # clone map, temporary directory, absolute path of input directory, and landmask
        self.cloneMap = iniItems.cloneMap
        self.tmpDir   = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']
        self.landmask = landmask

        attr = vos.getMapAttributesALL(self.cloneMap)
        self.nLat = int(attr['rows'])
        self.nLon = int(attr['cols'])

        # should we simulate hydrology during off-season
        self.OffSeason = bool(int(iniItems.globalOptions['OffSeason']))
        
        # soil parameters
        self.soil = parSoilAndTopo.SoilAndTopoParameters(iniItems, self.landmask)
        self.soil.read()
        # TODO: is there a better place to put this?
        if groundwater.WaterTable:
            self.soil.compute_capillary_rise_parameters()

        # crop parameters
        self.crop = parCrop.CropParameters(iniItems, self.landmask)
        self.crop.read()
        self.crop.computeCropCalendar(currTimeStep, meteo)
        
        # rotation parameters
        # - get nRotation
        self.rotationParameterFileNC = iniItems.rotationOptions['rotationParameterNC']
        self.read_rotation_parameters()

        # irrigation management parameters
        self.irrMgmtParameterFileNC = iniItems.irrMgmtOptions['irrMgmtParameterNC']
        self.read_irrigation_management_parameters()

        # field management parameters
        self.fieldMgmtParameterFileNC = iniItems.fieldMgmtOptions['fieldMgmtParameterNC']
        self.read_field_management_parameters()
        
        # TODO
        # # cellArea (unit: m2)
        # self.cellArea = vos.readPCRmapClone(iniItems.routingOptions['cellAreaMap'], \
        #                                     self.cloneMap,\
        #                                     self.tmpDir,\
        #                                     self.inputDir)
        # self.cellArea = pcr.ifthen(self.landmask, self.cellArea)

        # # TODO: rescale rotation fractions (necessary? - should sum to <= 1)

        # TODO: not sure whether this is the right place to initialize
        # DAP/SeasonCounter/GrowingSeason - are they re-initialized at each time step?
        

        # initialise season counter
        self.SeasonCounter = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))

        # initialise GDD counter

        # initialize TrRatio, rCor (? - what do these vars represent)
        self.TrRatio = np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        self.rCor = np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))

        """
        Counters
        """
        arr_zeros = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        self.AgeDays = arr_zeros
        self.AgeDays_NS = arr_zeros
        self.AerDays = arr_zeros
        self.IrrCum = arr_zeros
        self.DelayedGDDs = arr_zeros
        self.DelayedCDs = arr_zeros
        self.PctLagPhase = arr_zeros
        self.tEarlySen = arr_zeros
        self.GDDcum = arr_zeros
        self.DaySubmerged = np.zeros((self.nRotation, self.nLat, self.nLon))
        self.IrrNetCum = arr_zeros
        self.DAP = arr_zeros
        self.Epot = arr_zeros
        self.Tpot = arr_zeros
        self.GrowingSeason = arr_zeros

        """
        States (adapted from AOS_ReadModelInitialConditions lines 19-22)
        """
        self.PreAdj = np.full((self.crop.nCrop, self.nRotation, self.nLat, self.nLon), False)
        self.CropMature = np.full((self.crop.nCrop, self.nRotation, self.nLat, self.nLon), False)
        self.CropDead = np.full((self.crop.nCrop, self.nRotation, self.nLat, self.nLon), False)
        self.Germination = np.full((self.crop.nCrop, self.nRotation, self.nLat, self.nLon), False)
        self.PrematSenes = np.full((self.crop.nCrop, self.nRotation, self.nLat, self.nLon), False)
        self.HarvestFlag = np.full((self.crop.nCrop, self.nRotation, self.nLat, self.nLon), False)

        """
        Harvest index
        """
        arr_ones = np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        self.Stage = arr_ones
        self.Fpre = arr_ones
        self.Fpost = arr_ones
        self.fpost_dwn = arr_ones
        self.fpost_upp = arr_ones

        arr_zeros = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        self.HIcor_Asum = arr_zeros
        self.HIcor_Bsum = arr_zeros
        self.Fpol = arr_zeros
        self.sCor1 = arr_zeros
        self.sCor2 = arr_zeros

        """
        Crop growth initial conditions (adapted from AOS_ReadModelInitialConditions lines 39-55)
        """
        arr_zeros = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        self.CC = arr_zeros
        self.CCadj = arr_zeros
        self.CC_NS = arr_zeros
        self.CCadj_NS = arr_zeros
        self.Zroot = arr_zeros
        self.B = arr_zeros
        self.B_NS = arr_zeros
        self.HI = arr_zeros
        self.HIadj = arr_zeros
        self.CCxAct = arr_zeros
        self.CCxAct_NS = arr_zeros
        self.CCxW = arr_zeros
        self.CCxW_NS = arr_zeros
        self.CCxEarlySen = arr_zeros
        self.CCprev = arr_zeros
        self.CC0adj = arr_zeros  # TODO: see AOS_ReadModelInitialConditions - SeasonCounter???

        """
        Growth stage
        """
        self.GrowthStage = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))

        """
        Aeration stress (compartment level)
        """
        self.AerDaysComp = np.zeros((self.soil.nComp, self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        
        # Initialise variables for soil evaporation
        # TODO: may not need crop dimension
        self.Stage2 = np.full((self.nRotation, self.nLat, self.nLon), False)
        self.EvapZ = np.zeros((self.nRotation, self.nLat, self.nLon))
        self.Wstage2 = np.zeros((self.nRotation, self.nLat, self.nLon))
        self.Wsurf = np.zeros((self.nRotation, self.nLat, self.nLon))
        
        """
        Initial conditions
        """

        self.initialConditionFileNC = iniItems.landSurfaceOptions['initialConditionNC']

        # get the initial conditions (for every rotation)
        self.getInitialConditions(iniItems, initialState)
        
        # make iniItems available for the other methods/functions:
        self.iniItems = iniItems

    # def read_crop_parameters(self):

    #     # using LandCover.read_land_cover_parameters() as a template

    #     cropParams = ['CropType','CalendarType','crop.PlantingDate','crop.HarvestDate','SwitchGDD','Emergence','MaxRooting','Senescence','Maturity','HIstart','Flowering','YldForm','GDDmethod','PolHeatStress','PolColdStress','BioTempStress','PlantPop','Determinant','ETadj','LagAer','Tbase','Tupp','Tmax_up','Tmax_lo','Tmin_up','Tmin_lo','GDD_up','GDD_lo','fshape_b','PctZmin','Zmin','Zmax','fshape_r','fshape_ex','SxTopQ','SxBotQ','a_Tr','SeedSize','CCmin','CCx','CDC','CGC','Kcb','fage','WP','WPy','fsink','bsted','bface','HI0','HIini','dHI_pre','a_HI','b_HI','dHI0','exc','MaxFlowPct','p_up1','p_up2','p_up3','p_up4','p_lo1','p_lo2','p_lo3','p_lo4','fshape_w1','fshape_w2','fshape_w3','fshape_w4','Aer','beta','GermThr']

    #     for var in cropParams:
    #         vars(self)[var] = vos.netcdf2PCRobjCloneWithoutTime(self.cropParameterFileNC,\
    #                                                             var,\
    #                                                             cloneMapFileName=self.cloneMap)

    #     self.crop.nCrop = self.CropType.shape[0]

    def read_rotation_parameters(self):

        rotationParams = ['CropSequence']
        for var in rotationParams:
            vars(self)[var] = vos.netcdf2PCRobjCloneWithoutTime(self.rotationParameterFileNC,\
                                                                var,\
                                                                cloneMapFileName=self.cloneMap)

        self.nRotation = self.CropSequence.shape[1]
        
    def read_irrigation_management_parameters(self):

        irrMgmtParams = ['IrrMethod','IrrInterval','SMT1','SMT2','SMT3','SMT4','MaxIrr','AppEff','NetIrrSMT','WetSurf']
        for var in irrMgmtParams:
            vars(self)[var] = vos.netcdf2PCRobjCloneWithoutTime(self.irrMgmtParameterFileNC,\
                                                                var,\
                                                                cloneMapFileName=self.cloneMap)

        # create array of soil moisture target for the respective growth stages
        self.SMT = np.concatenate((self.SMT1[None,:], self.SMT2[None,:], self.SMT3[None,:], self.SMT4[None,:]), axis=0)

        # check if an irrigation schedule file is required
        if np.sum(self.IrrMethod == 3) > 0:
            if self.irrMgmtOptions['irrScheduleNC'] != "None":
                self.irrMgmtOptions['irrScheduleNC'] = vos.getFullPath(self.irrMgmtOptions[item], self.globalOptions['inputDir'])
                self.irrScheduleFileNC = irrMgmtOptions['irrScheduleNC']
            else:
                logger.error('IrrMethod equals 3 in some or all places, but irrScheduleNC is not set in configuration file')

        else:
            self.irrScheduleFileNC = None
        
    def read_field_management_parameters(self):

        fieldMgmtParams = ['Mulches','MulchPctGS','MulchPctOS','fMulch','Bunds','zBund','BundWater']
        for var in fieldMgmtParams:
            vars(self)[var] = vos.netcdf2PCRobjCloneWithoutTime(self.irrMgmtParameterFileNC,\
                                                                var,\
                                                                cloneMapFileName=self.cloneMap)            
            
    # TODO: make this a method of Crop object
    def getInitialConditions(self, iniItems, iniConditions = None):

        # TODO: see equivalent function in PCRGLOBWB - update rotation fractions

        # get initial conditions
        initialVars = ['th','Zroot']
        for var in initialVars:
            if iniConditions == None:
                # read input from NetCDF
                vars(self)[var] = vos.netcdf2PCRobjCloneWithoutTime(self.initialConditionFileNC,\
                                                                    var,\
                                                                    cloneMapFileName=self.cloneMap)
            else:
                vars(self)[var] = iniConditions[str(var)]

            # NB check this is doing what you think, because landmask is 2D
            # (lat,lon) and vars may be up to 4D (e.g. crop,rotation,lat,lon)
            vars(self)[var] = np.where(self.landmask, vars(self)[var], np.nan)

        # The following is adapted from AOS_ReadModelInitialConditions, lines 57-69:

        # "Get initial storage between surface bunds", if bunds are present.
        # Otherwise surface storage is zero.
        self.SurfaceStorage = np.zeros((self.nRotation, self.nLat, self.nLon))
        self.SurfaceStorageIni = np.zeros((self.nRotation, self.nLat, self.nLon))
        cond1 = ((self.Bunds == 0) & (self.zBund > 0.001))
        self.SurfaceStorage[cond1] = self.BundWater[cond1]
        cond2 = (cond1 & (self.SurfaceStorage > self.zBund))
        self.SurfaceStorage[cond1] = self.zBund[cond1]
        self.SurfaceStorageIni[cond1] = self.SurfaceStorage[cond1]
        
        # # - first, we set all aggregated states to zero (only the ones in mainStates): 
        # for var in self.mainStates: vars(self)[var] = pcr.scalar(0.0)
        # # - then we initiate them in the following loop of land cover types: 
        # for coverType in self.coverTypes:
        #     if iniConditions != None:
        #         self.landCoverObj[coverType].getICsLC(iniItems,iniConditions['landSurface'][coverType])
        #     else:
        #         self.landCoverObj[coverType].getICsLC(iniItems)
        #     # summarize/aggregate the initial states/storages (using the initial land cover fractions: previousFracVegCover)
        #     for var in self.mainStates:
        #         # - initial land cover fractions (dimensionless) 
        #         if isinstance(self.landCoverObj[coverType].previousFracVegCover, types.NoneType):
        #             self.landCoverObj[coverType].previousFracVegCover = self.landCoverObj[coverType].fracVegCover
        #         land_cover_fraction = self.landCoverObj[coverType].previousFracVegCover
        #         # - initial land cover states (unit: m)
        #         land_cover_states = vars(self.landCoverObj[coverType])[var]
        #         vars(self)[var]  += land_cover_states * land_cover_fraction

    def AOS_CheckGroundwaterTable(self, groundwater, currTimeStep):

        # expand soil properties to compartments
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_s  = self.soil.th_s[self.soil.layerIndex,:]
        
        # These are options which apply to the entire study area (not just a
        # selection of cells), so we can use an if statement.
        if groundwater.WaterTable:
        
            # get the mid point of each compartment (TODO: put this in initialization)
            zBot = np.cumsum(self.soil.dz)
            zTop = zBot - self.soil.dz
            zMid = (zTop + zBot) / 2
        
            # convenient to add compartment + rotation dimensions to groundwater
            # and rotation, latitude and longitude dimensions to zMid
            zGW = groundwater.zGW[None,None,:,:] * np.ones((self.soil.nComp,self.nRotation))[:,:,None,None]
            zMid = zMid[:,None,None,None] * np.ones((self.nRotation,self.nLat,self.nLon))[None,:,:,:]

            # "check if water table is within modelled soil profile"
            self.WTinSoilComp = (zMid >= zGW)
            self.th[self.WTinSoilComp] = th_s[self.WTinSoilComp]

            # Flatten WTinSoilComp to provide an array with dimensions
            # (nLat, nLon), indicating grid cells where the water table is in the
            # soil profile
            self.WTinSoil = np.sum(self.WTinSoilComp, axis=(0,1)).astype(bool)

            # "adjust compartment field capacity"
            self.soil.th_fc_adj = np.zeros((self.soil.nComp, self.nRotation, self.nLat, self.nLon))

            # get Xmax (TODO: find out what this variable represents)
            Xmax = np.zeros((self.soil.nComp,self.nRotation,self.nLat,self.nLon))
            cond1 = th_fc <= 0.1
            cond2 = th_fc >= 0.3
            cond3 = np.logical_not(cond1 | cond2) # i.e. 0.1 < fc < 0.3
            Xmax[cond1] = 1
            Xmax[cond2] = 2
            pF = 2 + 0.3 * (th_fc - 0.1) / 0.2
            Xmax_cond3 = np.exp(pF * np.log(10)) / 100
            Xmax[cond3] = Xmax_cond3[cond3]

            cond4 = (zGW < 0) | ((zGW - zMid) >= Xmax)

            # 'compartment' shows the index of the compartment to which each element
            # belongs (shallow -> deep, i.e. 1 is the shallowest)
            compartment = np.arange(1, self.soil.nComp + 1)[:,None,None,None] * np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:]

            # 'cond4_max_compartment' holds the index of the lowest compartment (i.e.
            # the maximum value) for which cond4 is met, cast to all compartments.
            # NB:
            # -> 'numpy.amax' returns the maximum value along an axis (0, or
            #    compartments, in this case)
            # -> multiply compartments by cond4 to set elements that do not equal the
            #    condition to zero, but retain the compartment number of elements that
            #    do meet the condition
            cond4_max_compartment = np.amax(compartment * cond4, axis=0)[None,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None]

            # now, identify compartments that are shallower than the deepest
            # compartment for which cond4 is met (array 'c'). NB:
            # -> because 'compartment' starts from 1, rather than zero, elements in
            #    'cond4_max_compartment' which are equal to zero (implying that none
            #    of the compartments meet cond4) are correctly set to False in array
            #    'cond4'
            cond4 = (compartment <= cond4_max_compartment)

            # 'cond4' is a special case because if ANY compartment meets the
            # condition then all overlying compartments are automatically assumed to
            # meet the condition. Thus in subsequent conditions we have to be careful
            # to ensure that True elements in 'cond4' do not also belong to 'cond5',
            # 'cond6' or 'cond7'. We use numpy.logical_not(...) for this purpose.
            cond5 = (th_fc >= th_s) & np.logical_not(cond4)
            cond6 = (zMid >= zGW) & np.logical_not(cond4 | cond5)
            cond7 = np.logical_not(cond4 | cond5 | cond6)
            dV = th_s - th_fc
            dFC = (dV / (Xmax ** 2)) * ((zMid - (zGW - Xmax)) ** 2)

            self.soil.th_fc_adj[cond4] = th_fc[cond4]
            self.soil.th_fc_adj[cond5] = th_fc[cond5]
            self.soil.th_fc_adj[cond6] = th_s[cond6]
            self.soil.th_fc_adj[cond7] = th_fc[cond7] + dFC[cond7]

        else:
            self.zGW = np.ones((self.nRotation, self.nLat, self.nLon)) * -999
            self.WTinSoil = np.full((self.nLat, self.nLon), False)
            self.soil.th_fc_adj = th_fc

    def PreIrrigation(self):

        # expand soil properties to compartments
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        
        # convenient to add crop dimension to th_wp, th_fc and current soil
        # moisture (th) - variables 'wp' and 'fc' have dimension (compartment, crop, rotation, lat, lon)
        th_wp = th_wp[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]
        th_fc = th_fc[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]
        th = self.th[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]
        
        # Zmin has dimension (crop,lat,lon), Zroot has dimension (crop,rotation,lat,lon)
        rootdepth = np.maximum(self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None], self.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100 # 2 decimal places

        # "determine critical water content threshold"
        # NetIrrSMT has dimension (crop, rotation, lat, lon)
        thCrit = (th_wp + (((self.NetIrrSMT[None,:,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None,None]) / 100) * (th_fc - th_wp)))

        # pre-irrigation requirement for each layer
        preIrr_req = (thCrit - th) * 1000 * self.soil.dz[:,None,None,None,None] * np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))[None,:,:,:,:]

        # calculate the various conditions under which pre-irrigation should be calculated
        cond1 = (self.IrrMethod == 4)[None,:,:,:,:] * np.full((self.soil.nComp), True)[:,None,None,None,None]
        cond2 = (self.DAP == 1)
        cond3 = ((rootdepth[None,:,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None,None]) >= (self.soil.dz[:,None,None,None,None] * np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))[None,:,:,:,:]))
        cond4 = (th < thCrit)

        # combine conditions
        cond5 = cond1 & cond2 & cond3 & cond4
        preIrr_req[np.logical_not(cond5)] = 0
        self.preIrr = np.sum(preIrr_req, axis=0) # compartment axis - dimensions are (crop, rotation, lat, lon)

    def Drainage(self):

        # expand soil properties to compartments
        ksat  = self.soil.ksat[self.soil.layerIndex,:]
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_s  = self.soil.th_s[self.soil.layerIndex,:]
        tau   = self.soil.tau[self.soil.layerIndex,:]
        
        # preallocate arrays
        drainsum = np.zeros((self.nRotation, self.nLat, self.nLon))
        thnew = np.zeros((self.soil.nComp, self.nRotation, self.nLat, self.nLon))
        self.FluxOut = np.zeros((self.soil.nComp, self.nRotation, self.nLat, self.nLon))

        for comp in range(self.soil.nComp):

            # ===========================
            # this section covers lines 16-50 of AOS_Drainage.m

            dthdt = np.zeros((self.nRotation, self.nLat, self.nLon))  # initialize

            # "Calculate drainage ability of compartment ii"
            cond1 = (self.th[comp,:] <= self.soil.th_fc_adj[comp,:])
            dthdt[cond1] = 0

            cond2 = (np.logical_not(cond1) & (self.th[comp,:] >= th_s[comp,:]))
            dthdt[cond2] = ((tau[comp,:] * th_s[comp,:]) - th_fc[comp,:])[cond2]

            cond3 = np.logical_not(cond1 | cond2)
            dthdt[cond3] = (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]) * ((np.exp(self.th[comp,:] - th_fc[comp,:]) - 1) / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond3]

            cond4 = ((cond2 | cond3) & ((self.th[comp,:] - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond4] = (self.th[comp,:] - self.soil.th_fc_adj[comp,:])[cond4]

            # "Drainage from compartment ii (mm)"
            draincomp = dthdt * self.soil.dz[comp]

            # "Check drainage ability of compartment ii against cumulative drainage from compartments above"
            excess = np.zeros((self.nRotation, self.nLat, self.nLon))
            prethick = self.soil.dzsum[comp] - self.soil.dz[comp]
            drainmax = dthdt * 1000 * prethick
            drainability = (drainsum <= drainmax)

            # ===========================
            # this section covers lines 52-65 of AOS_Drainage.m
            
            # "Drain compartment ii"
            # if drainability == True, "No storage needed. Update water content in compartment ii"
            thnew[comp,:][drainability] = (self.th[comp,:] - dthdt)[drainability]
            
            # "Update cumulative drainage (mm)"
            drainsum[drainability] += draincomp[drainability]

            # "Restrict cumulative drainage to saturated hydraulic conductivity and adjust excess drainage flow"
            cond6 = (drainability & (drainsum > ksat[comp,:]))
            excess[cond6] += (drainsum - ksat[comp,:])[cond6]
            drainsum[cond6] = ksat[comp,:][cond6]

            # ===========================
            # this section covers lines 66-84 of AOS_Drainage.m
            
            # if drainability == False, "Storage is needed"
            dthdt[np.logical_not(drainability)] = (drainsum / (1000 * prethick))[np.logical_not(drainability)]

            # "Calculate value of theta (thX) needed to provide a drainage ability equal to cumulative drainage
            thX = np.zeros((self.nRotation, self.nLat, self.nLon)) # preallocate
            cond7 = (np.logical_not(drainability) & (dthdt <= 0))
            thX[cond7] = self.soil.th_fc_adj[comp,:][cond7]

            cond8 = (np.logical_not(drainability) & np.logical_not(cond7) & (tau[comp,:] > 0))
            A = 1 + ((dthdt * (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)) / (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:])))
            thX[cond8] = (self.soil.th_fc_adj[comp,:] + np.log(A))[cond8]

            cond9 = (cond8 & (thX < self.soil.th_fc_adj[comp,:]))
            thX[cond9] = self.soil.th_fc_adj[comp,:][cond9]

            cond10 = (np.logical_not(drainability) & np.logical_not(cond7 | cond8))
            thX[cond10] = th_s[comp,:][cond10] + 0.01

            # ===========================
            # this section covers lines 85-122 of AOS_Drainage.m

            # "Check thX against hydraulic properties of current soil layer"
            cond11 = (np.logical_not(drainability) & (thX <= th_s[comp,:]))
            thnew[comp,:][cond11] = (self.th[comp,:] + (drainsum / (1000 * self.soil.dz[comp])))[cond11]

            # "Check updated water content against thX"
            cond12 = (cond11 & (thnew[comp,:] > thX))

            # "Cumulative drainage is the drainage difference between theta_x and new theta plus drainage ability at theta_x"
            drainsum[cond12] = ((thnew[comp,:] - thX) * 1000 * self.soil.dz[comp])[cond12]

            # "Calculate drainage ability for thX"
            cond13 = (cond12 & (thX <= self.soil.th_fc_adj[comp,:]))
            dthdt[cond13] = 0

            cond14 = (cond12 & np.logical_not(cond13) & (thX >= th_s[comp,:]))
            dthdt[cond14] = (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]))[cond14]

            cond15 = (cond14 & ((thX - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond15] = (thX - self.soil.th_fc_adj[comp,:])[cond15]

            cond16 = (cond12 & np.logical_not(cond13 | cond14))
            dthdt[cond16] = (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]) * ((np.exp(thX - th_fc[comp,:]) - 1) / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond16]

            cond17 = (cond16 & ((thX - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond17] = (thX - self.soil.th_fc_adj[comp,:])[cond17]

            # "Update drainage total"
            drainsum[cond12] += (dthdt * 1000 * self.soil.dz[comp])[cond12]

            # "Restrict cumulative drainage to saturated hydraulic conductivity and adjust excess drainage flow"
            cond18 = (cond12 & (drainsum > ksat[comp,:]))
            excess[cond18] += (drainsum - ksat[comp,:])[cond18]
            drainsum[cond18] = ksat[comp,:][cond18]

            # "Update water content"
            thnew[comp,:][cond12] = (thX - dthdt)[cond12]

            # ===========================
            # this section covers lines 123-151 of AOS_Drainage.m
            
            cond19 = (cond11 & np.logical_not(cond12) & (thnew[comp,:] > self.soil.th_fc_adj[comp,:]))

            # "Calculate drainage ability for updated water content"
            cond20 = (cond19 & (thnew[comp,:] <= self.soil.th_fc_adj[comp,:]))
            dthdt[cond20] = 0

            cond21 = (cond19 & np.logical_not(cond20) & (thnew[comp,:] >= th_s[comp,:]))
            dthdt[cond21] = (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]))[cond21]

            cond22 = (cond21 & ((thnew[comp,:] - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond22] = (thnew[comp,:] - self.soil.th_fc_adj[comp,:])[cond22]

            cond23 = (cond19 & np.logical_not(cond20 | cond21))
            dthdt[cond23] = (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]) * ((np.exp(thnew[comp,:] - th_fc[comp,:]) - 1) / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond23]

            cond24 = (cond23 & ((thnew[comp,:] - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond24] = (thnew[comp,:] - self.soil.th_fc_adj[comp,:])[cond24]

            # "Update water content in compartment ii"
            thnew[comp,:][cond19] = (thnew[comp,:] - dthdt)[cond19]

            # "Update cumulative drainage"
            drainsum[cond19] = (dthdt * 1000 * self.soil.dz[comp])[cond19]

            # "Restrict cumulative drainage to saturated hydraulic conductivity
            cond25 = (cond19 & (drainsum > ksat[comp,:]))
            excess[cond25] += (drainsum - ksat[comp,:])[cond25]
            drainsum[cond25] = ksat[comp,:][cond25]

            # ===========================
            # this section covers lines 152-156 of AOS_Drainage.m
            
            cond26 = (cond11 & np.logical_not(cond12 | cond19))

            # "Drainage and cumulative drainage are zero as water content has not risen above field capacity in compartment ii"
            drainsum[cond26] = 0

            # ===========================
            # this section covers lines 158-164 of AOS_Drainage.m

            cond27 = (np.logical_not(drainability) & np.logical_not(cond11)) # or, thX > th_s
            
            # "Increase water content in compartment ii with cumulative drainage from above"
            thnew[comp,:][cond27] = (self.th[comp,:] + (drainsum / (1000 * self.soil.dz[comp])))[cond27]

            # "Check new water content against hydraulic properties of soil layer"
            cond28 = (cond27 & (thnew[comp,:] <= th_s[comp,:]))

            # ===========================
            # this section covers lines 165-196 of AOS_Drainage.m
            
            cond29 = (cond28 & (thnew[comp,:] > self.soil.th_fc_adj[comp,:]))

            # "Calculate new drainage ability
            cond30 = (cond29 & (thnew[comp,:] <= th_fc[comp,:]))
            dthdt[cond30] = 0

            cond31 = (cond29 & np.logical_not(cond30) & (thnew[comp,:] >= th_s[comp,:]))
            dthdt[cond31] = (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]))[cond31]

            cond32 = (cond29 & np.logical_not(cond30 | cond31))
            dthdt[cond32] = (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]) * ((np.exp(thnew[comp,:] - th_fc[comp,:]) - 1) / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond32]

            cond33 = ((cond31 | cond32) & ((thnew[comp,:] - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond33] = (thnew[comp,:] - self.soil.th_fc_adj[comp,:])[cond33]

            # "Update water content in compartment ii"
            thnew[comp,:][cond29] -= (dthdt)[cond29]

            # "Update cumulative drainage"
            drainsum[cond29] = (dthdt * 1000 * self.soil.dz[comp])[cond29]
            
            # "Restrict cumulative drainage to saturated hydraulic conductivity and adjust excess drainage flow"
            cond35 = (cond29 & (drainsum > ksat[comp,:]))
            excess[cond35] += (drainsum - ksat[comp,:])[cond35]
            drainsum[cond35] = ksat[comp,:][cond35]

            drainsum[np.logical_not(cond29)] = 0

            # ===========================
            # this section covers lines 197-239 of AOS_Drainage.m

            cond36 = np.logical_not(cond28)

            # "Calculate excess drainage above saturation"
            excess[cond36] = ((thnew[comp,:] - th_s[comp,:]) * 1000 * self.soil.dz[comp])[cond36]

            # "Calculate drainage ability for updated water content"
            cond37 = (cond36 & (thnew[comp,:] <= self.soil.th_fc_adj[comp,:]))
            dthdt[cond37] = 0

            cond38 = (cond36 & (np.logical_not(cond37) & (thnew[comp,:] >= th_s[comp,:])))
            dthdt[cond38] = (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]))[cond38]

            cond39 = (cond36 & np.logical_not(cond37 | cond38))
            dthdt[cond39] = (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]) * ((np.exp(thnew[comp,:] - th_fc[comp,:]) - 1) / (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)))[cond39]

            cond40 = ((cond38 | cond39) & ((thnew[comp,:] - dthdt) < self.soil.th_fc_adj[comp,:]))
            dthdt[cond40] = (thnew[comp,:] - self.soil.th_fc_adj[comp,:])[cond40]

            # "Update water content in compartment ii"
            thnew[comp,:][cond36] = (th_s[comp,:] - dthdt)[cond36]

            # "Update drainage from compartment ii"
            draincomp[cond36] = (dthdt * 1000 * self.soil.dz[comp])[cond36]
        
            # "Update maximum drainage"
            drainmax[cond36] = (dthdt * 1000 * prethick)[cond36]

            # "Update excess drainage"
            cond42 = (cond36 & (drainmax > excess))
            drainmax[cond42] = excess[cond42]
            
            excess[cond36] -= drainmax[cond36]

            # "Update drainsum and restrict to saturated hydraulic conductivity of soil layer"
            drainsum[cond36] = (draincomp + drainmax)[cond36]

            cond43 = (cond36 & (drainsum > ksat[comp,:]))
            excess[cond43] += (drainsum - ksat[comp,:])[cond43]
            drainsum[cond43] = ksat[comp,:][cond43]

            # ===========================
            # this section covers lines 241-267 of AOS_Drainage.m

            # "Store output flux from compartment ii"
            self.FluxOut[comp,:] = drainsum

            # "Redistribute excess in compartment above"
            cond44 = (excess > 0)
            precomp = comp + 1
            while precomp != 0:
                
                # "Update compartment counter"
                precomp -= 1

                # "Update flux from compartment"
                if (precomp < comp):
                    self.FluxOut[precomp,:][cond44] -= excess[cond44]

                # "Increase water content to store excess"
                thnew[precomp,:][cond44] += (excess / (1000 * self.soil.dz[precomp]))[cond44]

                # "Limit water content to saturation and adjust excess counter"
                cond45 = (cond44 & (thnew[precomp,:] > th_s[precomp,:]))
                excess[cond45] = ((thnew[precomp,:] - th_s[precomp,:]) * 1000 * self.soil.dz[precomp])[cond45]
                thnew[precomp,:][cond45] = th_s[precomp,:][cond45]
                cond46 = (cond44 & np.logical_not(cond45))
                excess[cond46] = 0

        # "Update conditions and output"

        # "Total deep percolation (mm)
        self.DeepPerc = drainsum

        # "Water contents"
        self.th = thnew

    def RainfallPartition(self, meteo):

        # expand soil properties to compartments
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]
        
        # TODO: this seems to ignore interception - is this done elsewhere?
        
        # initialise runoff and infiltration arrays
        self.Runoff = np.zeros((self.nRotation,self.nLat,self.nLon))
        self.Infl = np.zeros((self.nRotation,self.nLat,self.nLon))

        # for convenience, cast dzsum to rotation, lat, lon
        dzsum = (self.soil.dzsum[:,None,None,None] * np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:])
        zcn = (self.soil.zCN[None,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None])

        cond1 = ((self.Bunds == 0) | self.zBund < 0.001)
        self.DaySubmerged[cond1] = 0

        cond2 = (cond1 & (self.soil.AdjCN == 1))

        # "Calulcate weighting factors by compartment"
        dzsum[dzsum > zcn] = zcn[dzsum > zcn]
        wx = 1.016 * (1 - np.exp(-4.16 * (dzsum / zcn)))

        # N.B. by the looks of things, in AquaCropOS xx is simply wx for the
        # overlying layer, with a value of zero for the top layer
        # N.B. wx[:-1,:,:,:] means all but the last element along the first dimension 
        xx = np.concatenate((np.zeros((self.nRotation, self.nLat, self.nLon))[None,:,:,:], wx[:-1,:,:,:]), axis=0)
        
        wrel = wx - xx
        wrel[wrel < 0] = 0
        wrel[wrel > 1] = 1

        # "Calculate relative wetness of topsoil"
        # TODO: should this update the main self.th???
        th = np.maximum(th_wp, self.th)

        wet_top_comp = wrel * ((th - th_wp) / (th_fc - th_wp))

        # cond3 replaces the comp_sto variable in AOS_RainfallPartition.m
        cond3 = dzsum >= zcn
        wet_top = np.sum(wet_top_comp * cond3, axis=0)
        wet_top[wet_top < 0] = 0
        wet_top[wet_top > 1] = 1

        cn = self.soil.CN
        cn[cond2] = np.round(self.soil.CNbot + (self.soil.CNtop - self.soil.CNbot) * wet_top)[cond2]

        # "Partition rainfall into runoff and infiltration"

        # for convenience
        P = meteo.precipitation[None,:,:] * np.ones((self.nRotation))[:,None,None]
        
        # first case - bunds present, so no runoff
        self.Runoff[np.logical_not(cond1)] = 0
        self.Infl[np.logical_not(cond1)] = P[np.logical_not(cond1)]

        S = (25400 / cn) - 254
        term = P - ((5 / 100) * S)
        
        cond4 = (cond1 & (term <= 0))
        self.Runoff[cond4] = 0
        self.Infl[cond4] = P[cond4]

        cond5 = (cond1 & np.logical_not(cond4))
        self.Runoff[cond5] = ((term ** 2) / (P + (1 - (5 / 100)) * S))[cond5]
        self.Infl[cond5] = (P - self.Runoff)[cond5]

    def RootZoneWater(self):
        # "Function to calculate actual and total available water in the root zone at current time step"

        # expand soil properties to compartments
        th_s   = self.soil.th_s[self.soil.layerIndex,:]
        th_fc  = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp  = self.soil.th_wp[self.soil.layerIndex,:]
        th_dry = self.soil.th_dry[self.soil.layerIndex,:]
        
        # convenient to add crop, rotation, lat, lon dimensions to dz and dzsum
        dz = (self.soil.dz[:,None,None,None,None] * np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))[None,:,:,:,:])
        dzsum = (self.soil.dzsum[:,None,None,None,None] * np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))[None,:,:,:,:])

        # "Calculate root zone water content and available water"
        # "Compartments covered by the root zone"

        # Zmin has dimension (crop,lat,lon), Zroot has dimension (crop,rotation,lat,lon)
        rootdepth = np.maximum(self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None], self.Zroot)
        rootdepth = np.round(rootdepth * 100) / 100

        # TODO: check (original Matlab code is "comp_sto = sum(Soil.Comp.dzsum < rootdepth) + 1")
        comp_sto = (dzsum < (rootdepth[None,:,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None,None]))

        # have to extend comp_sto by one (e.g. TTTFFF -> TTTTFF); easiest way is to
        # add an additional layer to the top and remove a layer from the bottom
        comp_sto_sum = np.sum(comp_sto, axis=0)
        comp_sto = np.concatenate((comp_sto_sum.astype(bool)[None,:,:,:,:], comp_sto[:-1,:]), axis=0)
        
        # "Fraction of compartment covered by root zone"
        # cond1 identifies compartments partially covered by root zone
        cond1 = (dzsum > ((rootdepth[None,:,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None,None]) / dz))

        factor = np.ones((self.soil.nComp, self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        factor[cond1] = (1 - ((dzsum - rootdepth[None,:,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None,None]) / dz))[cond1]
        factor[np.logical_not(comp_sto)] = 0

        # Now, 'factor' equals zero in all compartments NOT covered by the root
        # zone; thus Wr_comp, WrS_comp etc. will equal zero in the corresponding
        # compartments

        # "Actual water storage in root zone (mm)"
        Wr_comp = factor * 1000 * (self.th[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]) * dz
        Wr = np.sum(Wr_comp, axis=0)
        Wr[Wr < 0] = 0

        # "Water storage in root zone at saturation (mm)"
        WrS_comp = factor * 1000 * (th_s[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]) * dz
        WrS = np.sum(WrS_comp, axis=0)
        
        # "Water storage in root zone at field capacity (mm)"
        WrFC_comp = factor * 1000 * (th_fc[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]) * dz
        WrFC = np.sum(WrFC_comp, axis=0)

        # "Water storage in root zone at permanent wilting point (mm)"
        WrWP_comp = factor * 1000 * (th_wp[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]) * dz
        WrWP = np.sum(WrWP_comp, axis=0)
        
        # "Water storage in root zone at air dry (mm)"
        WrDry_comp = factor * 1000 * (th_dry[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]) * dz
        WrDry = np.sum(WrDry_comp, axis=0)

        # "Water storage in root zone at aeration stress threshold (mm)"
        WrAer_comp = factor * 1000 * ((th_s[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]) - (self.crop.Aer / 100)) * dz
        WrAer = np.sum(WrAer_comp, axis=0)

        # "Actual root zone water content (m3/m3)"
        self.thRZ_Act = Wr / (rootdepth * 1000)

        # "Root zone water content at saturation (m3/m3)"
        self.thRZ_Sat = WrS / (rootdepth * 1000)

        # "Root zone water content at field capacity (m3/m3)"
        self.thRZ_Fc  = WrFC / (rootdepth * 1000)

        # "Root zone water content at permanent wilting point (m3/m3)"
        self.thRZ_Wp  = WrWP / (rootdepth * 1000)

        # "Root zone water content at air dry (m3/m3)"
        self.thRZ_Dry = WrDry / (rootdepth * 1000)

        # "Root zone water content at aeration stress threshold (m3/m3)"
        self.thRZ_Aer = WrAer / (rootdepth * 1000)

        # "Calculate Total Available Water (mm)"
        self.TAW = WrFC - WrWP
        self.TAW[self.TAW < 0] = 0

        # "Calculate depletion (mm)"
        self.Dr = WrFC - Wr
        self.Dr[self.Dr < 0] = 0

    def Irrigation(self, meteo):

        # "Calculate root zone water content and depletion"
        self.RootZoneWater()
        
        # for convenience
        P = meteo.precipitation[None,:,:] * np.ones((self.nRotation))[:,None,None]

        # "Determine adjustment for inflows and outflows on current day"
        # dimensions of thRZ_* are (crop, rotation, lat, lon)
        cond1 = (self.thRZ_Act > self.thRZ_Fc)
        rootdepth = np.maximum(self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None], self.Zroot)
        AbvFc = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        AbvFc[cond1] = ((self.thRZ_Act - self.thRZ_Fc) * 1000 * rootdepth)[cond1]

        WCadj = self.Tpot + self.Epot - P + self.Runoff - AbvFc

        # "Determine irrigation depth (mm/day) to be applied

        # "Update growth stage if it is first day of a growing season
        cond1 = (self.GrowingSeason & (self.DAP == 1))
        self.GrowthStage[cond1] = 1

        # initialise Irr array
        self.Irr = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        
        # "Run irrigation depth calculation"

        # 1 - "Rainfed - no irrigation"
        cond2 = (self.GrowingSeason & (self.IrrMethod == 0))
        self.Irr[cond2] = 0

        # 2 - "Irrigation - soil moisture"
        cond3 = (self.GrowingSeason & (self.IrrMethod == 1))

        # "Get soil moisture target for current growth stage"
        growthstages = np.arange(1, 5)[:,None,None,None,None] * np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))[None,:,:,:,:]
        cond4 = (((self.GrowthStage[None,:,:,:,:] * np.ones((4))[:,None,None,None,None]) == growthstages) & self.GrowingSeason)
        SMT = np.zeros((4, self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        SMT[cond4] = self.SMT[cond4]
        SMT = np.sum(SMT, axis=0)

        # "Determine threshold to initiate irrigation"
        IrrThr = (1 - SMT / 100) * self.TAW
        
        # 3 - "Irrigation - fixed interval"
        cond5 = (self.GrowingSeason & (self.IrrMethod == 2))

        # "Get number of days in growing season so far (subtract 1 so that [we]
        # always irrigate first on day 1 of each growing season"
        nDays = self.DAP - 1

        # combine irrigation methods 1 and 2, because hereafter they are the same
        cond6 = (cond3 | cond5)
        
        # "Adjust depletion for inflows and outflows today"
        self.Dr[cond6] += WCadj[cond6]
        self.Dr[cond6] = np.maximum(self.Dr, 0)[cond6]

        # check if conditions for irrigation method 1 or 2 are met
        cond7 = (cond6 & ((self.Dr > IrrThr) | ((nDays % self.IrrInterval) == 0)))

        # if they are met, "irrigation occurs"
        IrrReq = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        IrrReq[cond7] = np.maximum(0, self.Dr)[cond7]

        # "Adjust irrigation requirements for application efficiency"
        EffAdj = ((100 - self.AppEff) + 100) / 100
        IrrReq[cond7] *= EffAdj[cond7]

        # "Limit irrigation to maximum depth"
        self.Irr[cond7] = np.minimum(self.MaxIrr, IrrReq)[cond7]

        # "Irrigation - pre-defined schedule"
        cond8 = (self.GrowingSeason & (self.IrrMethod == 3))

        # if irrScheduleFileNC is None then this method is not used
        if self.irrScheduleFileNC != None:
            IrrReq = vos.netcdf2PCRobjClone(self.irrScheduleFileNC,\
                                            "irrigation_schedule",\
                                            str(currTimeStep.fulldate),\
                                            useDoy = method_for_time_index,\
                                            cloneMapFileName = self.cloneMap,\
                                            LatitudeLongitude = True)
            self.Irr[cond8] = IrrReq[cond8]
        
        # "Irrigation - net irrigation" (included for clarity, but not strictly necessary)
        cond9 = (self.GrowingSeason & (self.IrrMethod == 4))

        # "Net irrigation calculation performed after transpiration, so
        # irrigation is zero here"
        self.Irr[cond9] = 0

        # "Update cumulative irrigation counter for growing season"
        self.IrrCum += self.Irr

    def Infiltration(self):
        # "Function to infiltrate incoming water"

        # expand soil properties to compartments
        ksat   = self.soil.ksat[self.soil.layerIndex,:]
        th_s   = self.soil.th_s[self.soil.layerIndex,:]
        th_fc  = self.soil.th_fc[self.soil.layerIndex,:]
        tau    = self.soil.tau[self.soil.layerIndex,:]
                
        # initialise variables
        ToStore = np.zeros((self.nRotation, self.nLat, self.nLon))
        RunoffIni = np.zeros((self.nRotation, self.nLat, self.nLon))
        
        # "Update infiltration rate for irrigation"
        # "Note: irrigation amount adjusted for specified application efficiency"

        # ***Line 11 of AOS_Infiltration.m***
        # sum along crop axis
        irr = self.Irr * (self.AppEff / 100)
        irr = np.sum(irr, axis=0)        
        self.Infl += irr

        # Bunds on field
        
        # Dimensions of Bunds is (rotation, lat, lon)
        cond1 = self.Bunds == 1  # bunds on field

        cond2 = (cond1 & (self.zBund > 0.001))
        InflTot = self.Infl + self.SurfaceStorage

        # ***Line 19 of AOS_Infiltration.m***
        # "Update surface storage and infiltration storage"
        cond3 = (cond2 & (InflTot > 0))
        cond4 = (cond3 & (InflTot > ksat[0,:]))

        # "Infiltration limited by saturated hydraulic conductivity of surface soil layer
        ToStore[cond4] = ksat[0,:][cond4]  # ksat of top soil layer

        # "Additional water ponds on surface
        self.SurfaceStorage[cond4] = (InflTot - ksat[0,:])[cond4]

        # ***Line 27 of AOS_Infiltration.m***
        # Otherwise all water infiltrates
        cond5 = (cond3 & np.logical_not(cond4))
        ToStore[cond5] = InflTot[cond5]

        # "Reset surface storage depth to zero
        self.SurfaceStorage[cond5] = 0
        
        # "Calculate additional runoff"
        cond6 = (cond3 & (self.SurfaceStorage > (self.zBund * 1000)))

        # "Water overtops bunds and runs off"
        RunoffIni[cond6] = (self.SurfaceStorage - (self.zBund * 1000))[cond6]

        # "Surface storage equal to bund height"
        self.SurfaceStorage[cond6] = (self.zBund * 1000)[cond6]

        # "No overtopping of bunds"
        cond7 = (cond3 & np.logical_not(cond6))
        RunoffIni[cond7] = 0
        
        cond8 = (cond2 & np.logical_not(cond3))

        # Bunds not on field

        # ***Line 49 of AOS_Infiltration.m***
        cond9 = (self.Bunds == 0)  # bunds not on field

        # "Infiltration limited by saturated hydraulic conductivity of top soil layer"
        cond10 = (cond9 & (self.Infl > ksat[0,:]))
        ToStore[cond10] = ksat[0,:][cond10]

        # "Additional water runs off"
        RunoffIni[cond10] = (self.Infl - ksat[0,:])[cond10]

        # Otherwise all water infiltrates
        cond11 = (cond9 & np.logical_not(cond10))
        ToStore[cond11] = self.Infl[cond11]
        RunoffIni[cond11] = 0   # included for clarity (RunoffIni already initialised as zero)

        # Initialise arrays
        thnew = np.zeros((self.soil.nComp, self.nRotation, self.nLat, self.nLon))
        
        # "Infiltrate incoming water"
        # ... TODO
        Runoff = np.zeros((self.nRotation, self.nLat, self.nLon))
        DeepPerc = np.zeros((self.nRotation, self.nLat, self.nLon))
        
        for comp in range(self.soil.nComp):

            # "Calculate saturated drainage ability"
            dthdtS = tau[comp,:] * (th_s[comp,:] - th_fc[comp,:])

            # "Calculate drainage factor"
            factor = ksat[comp,:] / (dthdtS * 1000 * self.soil.dz[comp])

            # "Calculate drainage ability required"
            dthdt0 = ToStore / (1000 * self.soil.dz[comp])  # TODO: are these the right units?

            theta0 = np.zeros((self.nRotation, self.nLat, self.nLon))  # initialise

            # "Check drainage ability
            cond12 = (dthdt0 < dthdtS)

            # "Calculate water content, thX, needed to meet drainage dthdt0"
            cond13 = (cond12 & (dthdt0 < dthdtS))
            theta0[cond13] = self.soil.th_fc_adj[comp,:][cond13]
            cond14 = (cond12 & np.logical_not(cond13))
            A = 1 + ((dthdt0 * (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)) / (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:])))
            theta0[cond14] = A[cond14]

            # "Limit thX to between saturation and field capacity"
            cond15 = (cond12 & (theta0 > th_s[comp,:]))
            theta0[cond15] = th_s[comp,:][cond15]
            cond16 = (cond12 & (theta0 < self.soil.th_fc_adj[comp,:]))
            theta0[cond16] = self.soil.th_fc_adj[comp,:][cond16]
            dthdt0[cond16] = 0

            # "Limit water content and drainage to saturation"
            cond17 = np.logical_not(cond12)
            theta0[cond17] = th_s[comp,:][cond17]
            dthdt0[cond17] = dthdtS[cond17]

            # "Calculate maximum water flow through compartment ii"
            drainmax = factor * dthdt0 * 1000 * self.soil.dz[comp]

            # "Calculate total drainage from compartment ii"
            drainage = drainmax + self.FluxOut[comp,:]

            # "Limit drainage to saturated hydraulic conductivity"
            cond18 = (drainage > ksat[comp,:])
            drainmax[cond18] = (ksat[comp,:] - self.FluxOut[comp,:])[cond18]

            # "Calculate difference between threshold and current water contents"
            diff = theta0 - self.th[comp,:]

            cond19 = (diff > 0)
            thnew[comp,:][cond19] = (thnew[comp,:] + (ToStore / (1000 * self.soil.dz[comp])))[cond19]
            cond20 = (cond19 & (thnew[comp,:] > theta0))

            # "Water remaining that can infiltrate to compartments below"
            ToStore[cond20] = ((thnew[comp,:] - theta0) * 1000 * self.soil.dz[comp])[cond20]
            thnew[comp,:][cond20] = theta0[cond20]

            # Otherwise, all infiltrating water has been stored
            cond21 = (cond19 & np.logical_not(cond20))
            ToStore[cond21] = 0

            # "Update outflow from current compartment (drainage + infiltration flows)
            self.FluxOut[comp,:] += ToStore

            # Calculate backup of water into compartments above
            excess = ToStore - drainmax
            excess[excess < 0] = 0

            # Update water to store
            ToStore -= excess

            # Redistribute excess to compartments above
            cond22 = (excess > 0)

            precomp = comp + 1
            while precomp != 0:

                # "Keep storing in compartments above until soil surface is reached"
                
                # "Update compartment counter"
                precomp -= 1

                # "Update outflow from compartment"
                self.FluxOut[precomp,:][cond22] = (self.FluxOut[precomp,:] - excess)[cond22]

                # "Update water content"
                thnew[precomp,:][cond22] += (excess / (1000 * self.soil.dz[precomp]))[cond22]

                # "Limit water content to saturation"
                cond23 = (cond22 & (thnew[precomp,:] > th_s[precomp,:]))
                # "Update excess to store"
                excess[cond23] = ((thnew[precomp,:] - th_s[precomp,:]) * 1000 * self.soil.dz[precomp])[cond23]
                thnew[precomp,:][cond23] = th_s[precomp,:][cond23]
                cond24 = (cond22 & np.logical_not(cond23))
                excess[cond24] = 0

            # Any leftover water not stored becomes runoff
            cond25 = (excess > 0)
            self.Runoff[cond25] += excess[cond25]
        
        # "Infiltration left to store after bottom compartment becomes deep
        # percolation (mm)
        self.DeepPerc = ToStore

        # TODO: check the above loop is updating DeepPerc and Runoff properly

        # "Update total runoff"
        Runoff += RunoffIni

        # "Update surface storage (if bunds are present)"
        cond26 = ((Runoff > RunoffIni) & (self.Bunds == 1) & self.zBund > 0.001)
        # "Increase surface storage"
        self.SurfaceStorage += (Runoff - RunoffIni)
        # "Limit surface storage to bund height"
        cond27 = (cond26 & (self.SurfaceStorage > (self.zBund * 1000)))
        # "Additional water above top of bunds becomes runoff"
        Runoff[cond27] = (RunoffIni + (self.SurfaceStorage - (self.zBund * 1000)))[cond27]
        # "Set surface storage to bund height"
        self.SurfaceStorage[cond27] = (self.zBund * 1000)[cond27]
        # Otherwise there is no additional overtopping of bunds
        cond28 = (cond26 & np.logical_not(cond27))
        Runoff[cond28] = RunoffIni[cond28]

        # "Store updated water contents"
        # TODO: check this has updated properly
        self.th = thnew

        # "Update deep percolation, surface runoff, and infiltration values"
        self.DeepPerc += DeepPerc
        self.Infl -= Runoff
        self.Runoff += Runoff

    def CapillaryRise(self, groundwater):

        ksat  = self.soil.ksat[self.soil.layerIndex,:]
        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]
        aCR   = self.soil.aCR[self.soil.layerIndex,:]
        bCR   = self.soil.bCR[self.soil.layerIndex,:]

        # "Get groundwater table elevation on current day"
        zGW = groundwater.zGW[None,:,:] * np.ones((self.nRotation))[:,None,None]

        # "Calculate capillary rise"
        if groundwater.WaterTable:

            # "Get maximum capillary rise for bottom compartment"
            zBot = np.sum(self.soil.dz)
            zBotMid = zBot - (self.soil.dz[-1] / 2)

            MaxCR = np.zeros((self.nRotation, self.nLat, self.nLon))  # initialise
            cond1 = ((ksat[-1,:] > 0) & (zGW > 0) & ((zGW - zBotMid) < 4))
            cond2 = (cond1 & (zBotMid >= zGW))
            cond3 = (cond1 & np.logical_not(cond2))
            MaxCR[cond2] = 99
            MaxCR[cond3] = (np.exp((np.log(zGW - zBotMid) - bCR[-1,:]) / aCR[-1,:]))[cond3]
            MaxCR[MaxCR > 99] = 99

            # "Find top of next soil layer that is not within modelled soil profile"
            # NB self.soil.layerIndex contains the index of the layer corresponding to
            # each compartment; hence the last value is the layer corresponding to the
            # bottom compartment.
            zTopLayer = np.cumsum(self.soil.zLayer[np.arange(0, self.soil.layerIndex[-1])])

            # "Check for restrictions on upward flow caused by properties of
            # compartments that are not modelled in the soil water balance"
            layeri = self.soil.layerIndex[-1] # index of layer corresponding to bottom compartment
            LimCR = np.zeros((self.nRotation, self.nLat, self.nLon))  # initialise
            while layeri < (self.soil.nLayer - 1):  # subtract 1 because of zero-based indexing
                layeri += 1                         # layer fully below bottom compartment
                cond4 = (zTopLayer < self.soil.zGW)
                cond5 = (cond4 & (self.soil.ksat[layeri,:] > 0) & (self.soil.zGW > 0) & ((self.soil.zGW - zTopLayer) < 4))
                cond6 = (cond5 & (zTopLayer >= self.soil.zGW))
                cond7 = (cond5 & np.logical_not(cond6))
                LimCR[cond6] = 99
                LimCR[cond7] = np.exp((np.log(self.soil.zGW - zTopLayer) - self.soil.bCR[layeri,:]) / self.soil.aCR[layeri,:])[cond7]
                LimCR[LimCR > 99] = 99
                
                cond8 = (cond4 & (MaxCR > LimCR))
                MaxCR[cond8] = LimCR[cond8]
                zTopLayer += self.soil.zLayer[layeri]

                
            # "Calculate capillary rise"
            compi = self.soil.nComp - 1  # "Start at bottom of root zone"
            WCr = np.zeros((self.nRotation, self.nLat, self.nLon))  # "Capillary rise counter"
            while compi >= 0:
                # "Proceed upwards until maximum capillary rise occurs, soil
                # surface is reached, or encounter a compartment where downward
                # drainage/infiltration has already occurred on current day"
                cond0 = ((np.round(MaxCR * 1000) > 0) & (np.round(self.FluxOut[compi,:] * 1000) == 0))
                Df = np.zeros((self.nRotation, self.nLat, self.nLon))
                cond1 = (cond0 & ((self.th[compi,:] >= th_wp[compi,:]) & (self.soil.fshape_cr > 0)))
                cond2 = (cond0 & np.logical_not(cond1))
                Df[cond1] = (1 - (((self.th[compi,:] - th_wp[compi,:]) / (self.soil.th_fc_adj[compi,:] - th_wp[compi,:])) ** self.soil.fshape_cr))[cond1]
                Df[Df > 1] = 1
                Df[Df < 0] = 0
                Df[cond2] = 1
                
                # "Calculate relative hydraulic conductivity"
                thThr = (th_wp[compi,:] + th_fc[compi,:]) / 2
                Krel = np.zeros((self.nRotation, self.nLat, self.nLon))
                cond3 = (cond0 & (self.th[compi,:] < thThr))
                cond4 = (cond3 & ((self.th[compi,:] < th_wp[compi,:]) | (thThr <= th_wp[compi,:])))
                cond5 = (cond3 & np.logical_not(cond4))
                Krel[cond5] = ((self.th[compi,:] - th_wp[compi,:]) / (thThr - th_wp[compi,:]))[cond5]
                cond6 = (cond0 & np.logical_not(cond3))
                Krel[cond6] = 1

                # "Check if room is available to store water from capillary rise"
                dth = np.zeros((self.nRotation, self.nLat, self.nLon))
                dth[cond0] = (self.soil.th_fc_adj[compi,:] - self.th[compi,:])[cond0]
                dth = np.round((dth * 1000) / 1000)

                # "Store water if room is available"
                dthMax = Krel * Df * MaxCR / (1000 * self.soil.dz[compi])
                CRcomp = np.zeros((self.nRotation, self.nLat, self.nLon))
                
                cond7 = (cond0 & ((dth > 0) & ((zBot - (self.soil.dz[-1] / 2)) < zGW)))
                cond8 = (cond7 & (dth >= dthMax))
                cond9 = (cond7 & np.logical_not(cond8))
                self.th[compi,:][cond8] += dthMax[cond8]
                CRcomp[cond8] = (dthMax * 1000 * self.soil.dz[compi])[cond8]
                MaxCR[cond8] = 0

                self.th[compi,:][cond9] = self.soil.th_fc_adj[compi,:][cond9]
                CRcomp[cond9] = (dth * 1000 * self.soil.dz[compi])[cond9]
                MaxCR[cond9] = ((Krel * MaxCR) - CRcomp)[cond9]

                WCr[cond7] += CRcomp[cond7]

                # "Update bottom elevation of compartment"
                zBot -= self.soil.dz[compi]

                # "Update compartment and layer counters"
                compi -= 1

                # "Update restriction on maximum capillary rise"
                if compi >= 0:
                    zBotMid = zBot - (self.soil.dz[compi] / 2)
                    cond10 = (cond0 & ((ksat[compi,:] > 0) & (zGW > 0) & ((zGW - zBotMid) < 4)))
                    cond11 = (cond10 & (zBotMid >= zGW))
                    cond12 = (cond10 & np.logical_not(cond11))
                    cond13 = (cond0 & np.logical_not(cond10))
                    LimCR[cond11] = 99
                    LimCR[cond12] = (np.exp((np.log(zGW - zBotMid) - bCR[compi,:]) / aCR[compi,:]))[cond12]
                    LimCR[LimCR > 99] = 99
                    LimCR[cond13] = 0
                    cond14 = (cond0 & (MaxCR > LimCR))
                    MaxCR[cond14] = LimCR[cond14]
            
            # "Store total depth of capillary rise
            self.CrTot = WCr

    def AOS_Germination(self):
        # "Function to check if crop has germinated"

        th_fc = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp = self.soil.th_wp[self.soil.layerIndex,:]

        dzsum = self.soil.dzsum[:,None,None,None,None] * np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))[None,:,:,:,:]
        dz = self.soil.dz[:,None,None,None,None] * np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))[None,:,:,:,:]
        zgerm = self.soil.zGerm[None,None,:,:,:] * np.ones((self.soil.nComp,self.crop.nCrop))[:,:,None,None,None]

        # Here we force zGerm to have a maximum value equal to the depth of the deepest soil compartment
        # TODO: confirm this is acceptable
        zgerm[zgerm > np.max(self.soil.dzsum)] = np.max(self.soil.dzsum)

        # "Find compartments covered by top soil layer affecting germination"
        # TODO: look at comp_sto in other submodules
        comp_sto = (dzsum >= zgerm)
        comp_idx = np.arange(1, self.soil.nComp + 1)[:,None,None,None,None] * np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))[None,:,:,:,:]
        comp_idx_rev = np.flip(comp_idx, axis=0)
        comp_idx_max = (self.soil.nComp + 1) - np.amax(comp_idx_rev * comp_sto, axis=0)
        comp_sto = (comp_idx <= (comp_idx_max[None,:,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None,None]))

        # "Calculate water content in top soil layer"
        Wr_comp   = np.zeros((self.soil.nComp, self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        WrFC_comp = np.zeros((self.soil.nComp, self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        WrWP_comp = np.zeros((self.soil.nComp, self.crop.nCrop, self.nRotation, self.nLat, self.nLon))

        # "Determine fraction of compartment covered by soil"
        factor = np.zeros((self.soil.nComp, self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        cond1 = (comp_sto & (dzsum > zgerm))
        cond2 = (comp_sto & np.logical_not(cond1))
        factor[cond1] = (1 - ((dzsum - zgerm) / dz))[cond1]
        factor[cond2] = 1
        factor *= (self.GrowingSeason[None,:,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None,None])

        # Now, 'factor' equals zero in all compartments NOT covered by the root
        # zone and/or not in growing season; thus Wr_comp, WrFC_comp, WrWP_comp
        # will be zero in those compartments

        # "Increment actual water storage (mm)
        Wr_comp = (factor * 1000 * (self.th[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]) * dz)
        Wr_comp[Wr_comp < 0] = 0
        Wr = np.sum(Wr_comp, axis=0)

        # "Increment water storage at field capacity (mm)"
        WrFC_comp = (factor * 1000 * (th_fc[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]) * dz)
        WrFC = np.sum(WrFC_comp, axis=0)

        # "Increment water storage at permanent wilting point (mm)"
        WrWP_comp = (factor * 1000 * (th_wp[:,None,:,:,:] * np.ones((self.crop.nCrop))[None,:,None,None,None]) * dz)
        WrWP = np.sum(WrWP_comp, axis=0)

        # "Calculate proportional water content"
        WcProp = 1 - ((WrFC - Wr) / (WrFC - WrWP))

        # "Check if water content is above germination threshold"
        germthr = (self.crop.GermThr[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])  # Germination threshold
        cond4 = ((WcProp >= germthr) & (np.logical_not(self.Germination)))
        self.Germination[cond4] = True

        # "Increment delayed growth time counters if germination is yet to occur"
        cond5 = (self.GrowingSeason & (np.logical_not(self.Germination)))
        self.DelayedCDs[cond5] += 1
        self.DelayedGDDs[cond5] += (self.crop.GDD[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[cond5]
                 
        # "Not in growing season so no germination calculation is performed"
        self.Germination[np.logical_not(self.GrowingSeason)] = False
        self.DelayedCDs[np.logical_not(self.GrowingSeason)] = 0
        self.DelayedGDDs[np.logical_not(self.GrowingSeason)] = 0

    def AOS_GrowthStage(self):
        # "Function to calculate number of growing degree days on current day"
        if self.crop.CalendarType == 1:
            tAdj = self.DAP - self.DelayedCDs
        elif self.crop.CalendarType == 2:
            tAdj = self.GDDcum - self.DelayedGDDs

        # self.GrowthStage = np.zeros((self.nCrop, self.nRotation, self.nLat, self.nLon))
            
        # "Update growth stage"
        cond1 = (self.GrowingSeason & (tAdj <= (self.crop.Canopy10Pct[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])))
        cond2 = (self.GrowingSeason & (tAdj <= (self.crop.MaxCanopy[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])))
        cond3 = (self.GrowingSeason & (tAdj <= (self.crop.Senescence[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])))
        cond4 = (self.GrowingSeason & (tAdj > (self.crop.Senescence[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])))

        self.GrowthStage[cond1] = 1
        self.GrowthStage[cond2] = 2
        self.GrowthStage[cond3] = 3
        self.GrowthStage[cond4] = 4

    def AOS_RootDevelopment(self, groundwater):
        # "Function to calculate root zone expansion"

        # "If today is the first day of season, root depth is equal to minimum depth"
        cond1 = (self.GrowingSeason & (self.DAP == 1))
        self.Zroot[cond1] = (self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[cond1]

        # "Adjust time for any delayed development"
        if self.crop.CalendarType == 1:
            tAdj = self.DAP - self.DelayedCDs
        elif self.crop.CalendarType == 2:
            tAdj = self.GDDcum - self.DelayedGDDs

        # "Calculate root expansion"
        Zini = self.crop.Zmin * (self.crop.PctZmin / 100)
        t0 = np.round(self.crop.Emergence / 2)
        tmax = self.crop.MaxRooting

        if self.crop.CalendarType == 1:
            tOld = tAdj - 1
        elif self.crop.CalendarType == 2:
            tOld = tAdj - (self.crop.GDD[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])

        # "Potential root depth on previous day"
        ZrOld = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        cond2 = (self.GrowingSeason & (tOld >= tmax))
        cond3 = (self.GrowingSeason & (tOld <= t0))
        cond4 = (self.GrowingSeason & (np.logical_not(cond2 | cond3)))
        ZrOld[cond2] = (self.crop.Zmax[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[cond2]
        ZrOld[cond3] = (Zini[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[cond3]
        X = (tOld - t0) / (tmax - t0)
        exponent = 1 / (self.crop.fshape_r[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])  # exponent
        ZrOld[cond4] = (((Zini + (self.crop.Zmax - Zini))[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]) * np.power(X, exponent))[cond4]
        cond5 = (self.GrowingSeason & (ZrOld < (self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])))
        ZrOld[cond5] = (self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[cond5]

        # Potential root depth on current day
        Zr = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        cond6 = (self.GrowingSeason & (tAdj >= tmax))
        cond7 = (self.GrowingSeason & (tAdj <= t0))
        cond8 = (self.GrowingSeason & (np.logical_not(cond6 | cond7)))
        Zr[cond6] = (self.crop.Zmax[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[cond6]
        Zr[cond7] = (Zini[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[cond7]
        X = (tAdj - t0) / (tmax - t0)
        Zr[cond8] = (((Zini + (self.crop.Zmax - Zini))[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]) * np.power(X, exponent))[cond8]

        cond9 = (self.GrowingSeason & (Zr < (self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])))
        Zr[cond9] = (self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[cond9]

        # "Determine rate of change"
        dZr = Zr - ZrOld

        # "Adjust rate of expansion for any stomatal water stress"
        cond10 = (self.GrowingSeason & (self.TrRatio < 0.9999))
        cond101 = (cond10 & (self.crop.fshape_ex[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None] >= 0))
        cond102 = (cond10 & np.logical_not(cond101))
        dZr[cond101] = (dZr * self.TrRatio)[cond101]
        fAdj = (np.exp(self.TrRatio * self.crop.fshape_ex[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]) - 1) / (np.exp(self.crop.fshape_ex[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]) - 1)
        dZr[cond102] = (dZr * fAdj)[cond102]

        # "Adjust root expansion for failure to germinate (roots cannot expand
        # if crop has not germinated)"
        dZr[np.logical_not(self.Germination)] = 0

        # "Get new rooting depth"
        self.Zroot[self.GrowingSeason] = (self.Zroot + dZr)[self.GrowingSeason]

        # "Adjust root depth if restrictive soil layer is present that limits depth of root expansion"
        zres = (self.soil.zRes[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None])
        sxtop = (self.crop.SxTop[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])
        sxbot = (self.crop.SxBot[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])
        
        cond11 = (self.GrowingSeason & (zres > 0))
        cond111 = (cond11 & (self.Zroot > zres))
        self.rCor[cond111] = ((2 * (self.Zroot / zres) * ((sxtop + sxbot) / 2) - sxtop) / sxbot)[cond111]
        self.Zroot[cond111] = zres[cond111]

        # "Limit rooting depth if groundwater table is present (roots cannot develop below the water table)"
        zGW = groundwater.zGW[None,None,:,:] * np.ones((self.crop.nCrop, self.nRotation))[:,:,None,None]
        cond12 = (self.GrowingSeason & groundwater.WaterTable & (zGW > 0))
        cond121 = (cond12 & (self.Zroot > zGW))
        self.Zroot[cond121] = zGW[cond121]
        cond122 = (cond121 & (self.Zroot < self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]))
        self.Zroot[cond122] = (self.crop.Zmin[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[cond122]

        self.Zroot[np.logical_not(self.GrowingSeason)] = 0

    def AOS_WaterStress(self, meteo, beta):
        # "Function to calculate water stress coefficients"

        # "Calculate relative root zone water depletion for each stress type"
        nstress = 4
        p_up = np.concatenate((self.crop.p_up1[None,:],\
                               self.crop.p_up2[None,:],\
                               self.crop.p_up3[None,:],\
                               self.crop.p_up4[None,:]), axis=0)

        p_lo = np.concatenate((self.crop.p_lo1[None,:],\
                               self.crop.p_lo2[None,:],\
                               self.crop.p_lo3[None,:],\
                               self.crop.p_lo4[None,:]), axis=0)

        fshape_w = np.concatenate((self.crop.fshape_w1[None,:],\
                                   self.crop.fshape_w2[None,:],\
                                   self.crop.fshape_w3[None,:],\
                                   self.crop.fshape_w4[None,:]), axis=0)
        
        p_up = p_up[:,:,None,:,:] * np.ones((self.nRotation))[None,None,:,None,None]
        p_lo = p_lo[:,:,None,:,:] * np.ones((self.nRotation))[None,None,:,None,None]
        fshape_w = fshape_w[:,:,None,:,:] * np.ones((self.nRotation))[None,None,:,None,None]
        
        # "Adjust stress thresholds for Et0 on current day (don't do this for
        # pollination water stress coefficient)"

        # TODO: in the future, ET0 may be provided for each rotation, by downscaling temperature
        
        et0 = meteo.referencePotET[None,None,:,:] * np.ones((self.crop.nCrop,self.nRotation))[:,:,None,None]
        cond1 = (self.crop.ETadj == 1)[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        cond1 = cond1.astype(bool)
        for stress in range(3):
            p_up[stress,:][cond1] = (p_up[stress,:] + (0.04 * (5 - et0)) * (np.log10(10 - 9 * p_up[stress,:])))[cond1]
            p_lo[stress,:][cond1] = (p_lo[stress,:] + (0.04 * (5 - et0)) * (np.log10(10 - 9 * p_lo[stress,:])))[cond1]

        # "Adjust senescence threshold if early senescence triggered"
        if beta:
            cond2 = (self.tEarlySen > 0)  # dimensions of tEarlySen are crop,rotation,lat,lon
            beta = (self.crop.beta[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])
            p_up[2,:][cond2] = (p_up[2,:] * (1 - (beta / 100)))[cond2]

        # "Limit values"
        p_up[p_up < 0] = 0
        p_lo[p_lo < 0] = 0
        p_up[p_up > 1] = 1
        p_lo[p_lo > 1] = 1

        # "Calculate relative depletion"
        Drel = np.zeros((nstress, self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        cond1 = (self.Dr <= (p_up * self.TAW))
        cond2 = (self.Dr >  (p_up * self.TAW)) & (self.Dr < (p_lo * self.TAW))
        cond3 = (self.Dr >= (p_lo * self.TAW))
        Drel[cond1] = 0         # no water stress
        Drel[cond2] = (1 - ((p_lo - (self.Dr / self.TAW)) / (p_lo - p_up)))[cond2]  # partial water stress
        Drel[cond3] = 1         # full water stress

        # "Calculate root zone stress coefficients"
        Ks = (1 - ((np.exp(Drel[np.arange(0,3),:] * fshape_w[np.arange(0,3),:]) - 1) / (np.exp(fshape_w[np.arange(0,3),:]) - 1)))
        # "Water stress coefficient for leaf expansion"
        self.Ksw_Exp = Ks[0,:]
        # "Water stress coefficient for stomatal closure"
        self.Ksw_Sto = Ks[1,:]
        # "Water stress coefficient for senescence"
        self.Ksw_Sen = Ks[2,:]
        # "Water stress coefficient for pollination failure"
        self.Ksw_Pol = 1 - Drel[3,:]
        # "Mean water stress coefficient for stomatal closure"
        self.Ksw_StoLin = 1 - Drel[1,:]

    def AOS_CCDevelopment(self, CC0, CCx, CGC, CDC, dt, Mode):
        # "Function to calculate canopy cover development by end of the current simulation day"
        CC = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        if Mode == 'Growth':
            # "Exponential growth stage"
            CC[self.GrowingSeason] = (CC0 * np.exp(CGC * dt))[self.GrowingSeason]
            cond1 = (self.GrowingSeason & (CC > (CCx / 2)))
            # "Exponential decay stage"
            CC[cond1] = (CCx - 0.25 * (CCx / CC0) * CCx * np.exp(-CGC * dt))[cond1]
            # "Limit CC to CCx"
            cond2 = (self.GrowingSeason & (CC > CCx))
            CC[cond2] = CCx[cond2]
        elif Mode == "Decline":
            cond3 = (self.GrowingSeason & (CCx < 0.001))
            cond4 = (self.GrowingSeason & np.logical_not(cond3))
            CC[cond3] = 0
            CC[cond4] = (CCx * (1 - 0.05 * (np.exp(dt * (CDC / CCx)) - 1)))[cond4]

        # "Limit canopy cover to between 0 and 1
        CC = np.clip(CC, 0, 1)
        return CC

    def AOS_CCRequiredTime(self, CCprev, CC0, CCx, CGC, CDC, dt, tSum, Mode):
        # "Function to find required time to reach CC at end of previous day, given current CGC or CDC"

        # "Get CGC and/or time (GDD or CD) required to reach CC on previous day"
        if Mode == 'CGC':
            CGCx = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
            cond1 = (CCprev <= (CCx / 2))
            cond2 = np.logical_not(cond1)
            CGCx[cond1] = ((np.log(CCprev / CC0)) / (tSum - dt))[cond1]
            CGCx[cond2] = ((np.log((0.25 * CCx * CCx / CC0) / (CCx - CCprev))) / (tSum - dt))[cond2]
            tReq = (tSum - dt) * (CGCx / CGC)
        elif Mode == 'CDC':
            tReq = ((np.log(1 + (1 - CCprev / CCx) / 0.05)) / (CDC / CCx))

        return tReq
                    
    def AOS_AdjustCCx(self, CCprev, CC0, CCx, CGC, CDC, dt, tSum):
        # "Function to adjust CCx value for changes in CGC due to water stress during the growing season"

        # "Get time required to reach CC on previous day"
        tCCtmp = self.AOS_CCRequiredTime(CCprev, CC0, CCx, CGC, CDC, dt, tSum, 'CGC')

        CanopyDevEnd = self.crop.CanopyDevEnd[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        
        # "Determine CCx adjusted"
        cond1 = (tCCtmp > 0)
        tCCtmp[cond1] += ((CanopyDevEnd - tSum) + dt)[cond1]
        CCxAdj_tmp = self.AOS_CCDevelopment(CC0, CCx, CGC, CDC, tCCtmp, 'Growth')
        CCxAdj = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        CCxAdj[cond1] = CCxAdj_tmp[cond1]
        return CCxAdj

    def AOS_UpdateCCxCDC(self, CCprev, CDC, CCx, dt):
        # "Function to update CCx and CDC parameter values for rewatering in late
        # season of an early declining canopy"

        # "Get adjusted CCx"
        CCxAdj = CCprev / (1 - 0.05 * (np.exp(dt * (CDC / CCx)) - 1))

        # "Get adjusted CDC"
        CDCadj = CDC * (CCxAdj / CCx)

        return CCxAdj,CDCadj

    def AOS_CanopyCover(self, meteo):
        """
        Function to simulate canopy growth/decline
        """

        # TODO: Check against NewCond/InitCond in original Matlab code

        # "Calculate root zone water content"
        self.RootZoneWater()

        # "Determine if water stress is occurring"
        self.AOS_WaterStress(meteo, beta = True)

        # Add rotation to dimensions of crop parameters (i.e. crop/lat/lon -> crop/rotation/lat/lon)
        Emergence = self.crop.Emergence[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        Maturity = self.crop.Maturity[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        CanopyDevEnd = self.crop.CanopyDevEnd[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        CC0 = self.crop.CC0[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        CCx = self.crop.CCx[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        CDC = self.crop.CDC[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        CGC = self.crop.CGC[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        Senescence = self.crop.Senescence[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        GDD = self.crop.GDD[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]

        # Preallocate some variables
        CCxAdj = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        CDCadj = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        CCsen = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))

        # Save initial values of certain variables
        self.CCprev = self.CC
        
        """
        "Get canopy cover growth over time"
        """
        if self.crop.CalendarType == 1:
            tCC = self.DAP
            dtCC = np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
            tCCadj = self.DAP - self.DelayedCDs
        elif self.crop.CalendarType == 2:
            tCC = self.GDDcum
            dtCC = GDD
            tCCadj = self.GDDcum - self.DelayedGDDs

        """
        "Canopy development (potential)"
        """
        # "No canopy development before emergence/germination or after maturity"
        cond1 = (self.GrowingSeason & ((tCC < Emergence) | (np.round(Emergence) > Maturity)))
        self.CC_NS[cond1] = 0

        # "Canopy growth can occur"
        cond2 = (self.GrowingSeason & (tCC < CanopyDevEnd))

        # "Very small initial CC as it is first day or due to senescence. In this
        # case assume no leaf expansion stress"
        cond21 = (cond2 & (self.CC_NS <= CC0))
        self.CC_NS[cond21] = (CC0 * np.exp(CGC * dtCC))[cond21]
        
        # "Canopy growing"
        cond22 = (cond2 & np.logical_not(cond21))
        tmp_tCC = tCC - Emergence
        tmp_CC = self.AOS_CCDevelopment(CC0, CCx, CGC, CDC, tmp_tCC, 'Growth')
        self.CC[cond22] = tmp_CC[cond22]

        # "Update maximum canopy cover size in growing season"
        self.CCxAct_NS[cond2] = self.CC_NS[cond2]
        
        # "No more canopy growth is possible or canopy in decline"
        cond3 = (self.GrowingSeason & (tCC > CanopyDevEnd))
        # "Set CCx for calculation of withered canopy effects"
        self.CCxW_NS[cond3] = self.CCxAct_NS[cond3]
        
        # "Mid-season stage - no canopy growth"
        cond31 = (cond3 & (tCC < Senescence))
        # TODO: original Matlab code is NewCond.CC_NS = InitCond.CC_NS
        self.CCxAct_NS[cond31] = self.CC_NS[cond31]

        # "Late-season stage - canopy decline"
        cond32 = (cond3 & np.logical_not(cond31))
        tmp_tCC = tCC - Emergence
        tmp_CC = self.AOS_CCDevelopment(CC0, CCx, CGC, CDC, tmp_tCC, 'Decline')
        self.CC[cond32] = tmp_CC[cond32]

        """
        "Canopy development (actual)
        """
        # "No canopy development before emergence/germination or after maturity"
        cond4 = (self.GrowingSeason & ((tCCadj < Emergence) | (np.round(tCCadj) > Maturity)))
        self.CC[cond4] = 0

        # "Canopy growth can occur"
        cond5 = (self.GrowingSeason & (tCCadj < CanopyDevEnd) & np.logical_not(cond4))
        cond51 = (cond5 & (self.CCprev <= self.CC0adj))

        # "Very small initial CC as it is first day or due to senescence. In
        # this case, assume no leaf expansion stress"
        self.CC[cond51] = (self.CC0adj * np.exp(CGC * dtCC))[cond51]

        # "Canopy growing"        
        cond52 = (cond5 & np.logical_not(cond51))

        # "Canopy approaching maximum size"
        cond521 = (cond52 & (self.CCprev >= (0.9799 * CCx)))
        tmp_tCC = tCC - Emergence
        tmp_CC = self.AOS_CCDevelopment(CC0, CCx, CGC, CDC, tmp_tCC, 'Growth')
        self.CC[cond521] = tmp_CC[cond521]
        self.CC0adj[cond521] = CC0[cond521]

        # ***Line 83 of AOS_CanopyCover.m***
        
        # "Adjust canopy growth coefficient for leaf expansion water stress effects"
        cond522 = (cond52 & np.logical_not(cond521))
        CGCadj = CGC * self.Ksw_Exp

        # "Adjust CCx for change in CGC
        cond5221 = (cond522 & (CGCadj > 0))
        tmp_CCxAdj = self.AOS_AdjustCCx(self.CCprev, self.CC0adj, CCx, CGCadj, CDC, dtCC, tCCadj)
        CCxAdj[cond5221] = tmp_CCxAdj[cond5221]

        cond52211 = (cond5221 & (CCxAdj > 0))

        # "Approaching maximum canopy size"
        cond522111 = (cond52211 & (np.abs(self.CCprev - CCx) < 0.00001))
        tmp_tCC = tCC - Emergence
        tmp_CC = self.AOS_CCDevelopment(CC0, CCx, CGC, CDC, tmp_tCC, 'Growth')
        self.CC[cond522111] = tmp_CC[cond522111]

        # "Determine time required to reach CC on previous day, given CGCadj value"
        cond522112 = (cond52211 & np.logical_not(cond522111))
        tReq = self.AOS_CCRequiredTime(self.CCprev, self.CC0adj, CCxAdj, CGCadj, CDC, dtCC, tCCadj, 'CGC')

        # "Calculate GDD's for canopy growth"
        tmp_tCC = tReq + dtCC

        # "Determine new canopy size"
        cond5221121 = (cond522112 & (tmp_tCC > 0))
        tmp_CC = self.AOS_CCDevelopment(self.CC0adj, CCxAdj, CGCadj, CDC, tmp_tCC, 'Growth')
        self.CC[cond5221121] = tmp_CC[cond5221121]

        # ***Line 118 of AOS_CanopyCover.m***
        
        # "No canopy growth"
        cond5222 = (cond522 & np.logical_not(cond5221))

        # "Update CC0 if current canopy cover if less than initial canopy cover
        # size at planting"
        cond52221 = (cond5222 & (self.CC < self.CC0adj))
        self.CC0adj[cond52221] = self.CC[cond52221]

        # "Update actual maximum canopy cover size during growing season"
        cond53 = (cond5 & (self.CC > self.CCxAct))
        self.CCxAct[cond53] = self.CC[cond53]

        # ***Line 132 of AOS_CanopyCover.m***
        
        # "No more canopy growth is possible or canopy is in decline"
        cond6 = (self.GrowingSeason & (tCCadj > CanopyDevEnd))

        # "Mid-season stage - no canopy growth"
        cond61 = (cond6 & (tCCadj < Senescence))
        self.CC[cond61] = self.CCprev[cond61]

        # "Update actual maximum canopy cover size during growing season"
        cond611 = (cond61 & (self.CC > self.CCxAct))
        self.CCxAct[cond611] = self.CC[cond611]

        # "Late season stage - canopy decline"
        # "Adjust canopy decline coefficient for difference between actual and potential CCx"
        cond62 = (cond6 & np.logical_not(cond61))
        CDCadj[cond62] = (CDC * (self.CCxAct / CCx))[cond62]

        # "Determine new canopy size"
        tmp_tCC = tCCadj - Senescence
        tmp_CC = self.AOS_CCDevelopment(self.CC0adj, self.CCxAct, CGC, CDCadj, tmp_tCC, 'Decline')
        self.CC[cond62] = tmp_CC[cond62]

        # "Check for crop growth termination"
        # "Crop has died"
        cond63 = (cond6 & ((self.CC < 0.001) & np.logical_not(self.CropDead)))
        self.CC[cond63] = 0
        self.CropDead[cond63] = True

        # ***Line 160 of AOS_CanopyCover.m***
        
        """
        Canopy senescence due to water stress (actual)
        """
        cond7 = (self.GrowingSeason & (tCCadj >= Emergence))

        # "Check for early canopy senescence starting/continuing due to severe
        # water stress"
        cond71 = (cond7 & ((tCCadj < Senescence) | (self.tEarlySen > 0)))

        # "Early canopy senescence"
        cond711 = (cond71 & (self.Ksw_Sen < 1))
        self.PrematSenes[cond711] = True
        # "No prior early senescence"
        cond7111 = (cond711 & (self.tEarlySen == 0))
        self.CCxEarlySen[cond7111] = self.CCprev[cond7111]

        # "Increment early senescence GDD counter"
        self.tEarlySen[cond711] += dtCC[cond711]

        # "Adjust canopy decline coefficient for water stress"
        self.AOS_WaterStress(meteo, beta = False)
        cond7112 = (cond711 & (self.Ksw_Sen > 0.99999))
        CDCadj[cond7112] = 0.0001
        cond7113 = (cond711 & np.logical_not(cond7112))
        CDCadj[cond7113] = ((1 - (self.Ksw_Sen ** 8)) * CDC)[cond7113]

        # "Get new canopy cover size after senescence"
        cond7114 = (cond711 & (self.CCxEarlySen < 0.001))
        CCsen[cond7114] = 0
        cond7115 = (cond711 & np.logical_not(cond7114))

        # "Get time required to reach CC at end of previous day, given CDCadj"
        tReq = self.AOS_CCRequiredTime(self.CCprev, self.CC0adj, self.CCxEarlySen, CGC, CDCadj, dtCC, tCCadj, 'CDC')

        # "Calculate GDD's for canopy decline"
        tmp_tCC = tReq + dtCC

        # "Determine new canopy size"
        tmp_CCsen = self.AOS_CCDevelopment(self.CC0adj, self.CCxEarlySen, CGC, CDCadj, tmp_tCC, 'Decline')
        CCsen[cond7115] = tmp_CCsen[cond7115]

        # "Update canopy cover size"
        cond7116 = (cond711 & (tCCadj < Senescence))

        # "Limit CC to CCx"
        CCsen[cond7116] = np.clip(CCsen, None, CCx)[cond7116]

        # "CC cannot be greater than value on previous day"
        self.CC[cond7116] = np.clip(self.CC, None, self.CCprev)[cond7116]

        # "Update maximum canopy cover size during growing season"
        self.CCxAct[cond7116] = self.CC[cond7116]

        # "Update CC0 if current CC is less than initial canopy cover size at planting"
        cond71161 = (cond7116 & (self.CC < CC0))
        self.CC0adj[cond71161] = self.CC[cond71161]
        cond71162 = (cond7116 & np.logical_not(cond71161))
        self.CC0adj[cond71162] = CC0[cond71162]

        # "Update CC to account for canopy cover senescence due to water stress"
        cond7117 = (cond711 & np.logical_not(cond7116))
        cond71171 = (cond7117 & (CCsen < self.CC))
        self.CC[cond71171] = CCsen[cond71171]

        # "Check for crop growth termination
        cond7118 = (cond711 & ((self.CC < 0.001) & np.logical_not(self.CropDead)))
        self.CC[cond7118] = 0
        self.CropDead[cond7118] = True

        # "No water stress"
        cond712 = (cond71 & np.logical_not(cond711))
        self.PrematSenes[cond712] = False

        # "Rewatering of canopy in late season"
        cond7121 = (cond712 & ((tCCadj > Senescence) & (self.tEarlySen > 0)))

        # "Get new values for CCx and CDC"
        tmp_tCC = tCCadj - dtCC - Senescence
        tmp_CCxAdj,tmp_CDCadj = self.AOS_UpdateCCxCDC(self.CCprev, CDC, CCx, tmp_tCC)
        CCxAdj[cond7121] = tmp_CCxAdj[cond7121]
        CDCadj[cond7121] = tmp_CDCadj[cond7121]

        # "Get new CC value for end of current day
        tmp_tCC = tCCadj - Senescence
        tmp_CC = self.AOS_CCDevelopment(self.CC0adj, CCxAdj, CGC, CDCadj, tmp_tCC, 'Decline')
        self.CC[cond7121] = tmp_CC[cond7121]

        # "Check for crop growth termination"
        cond71211 = (cond7121 & ((self.CC < 0.001) & np.logical_not(self.CropDead)))
        self.CC[cond71211] = 0
        self.CropDead[cond71211] = True

        # "Reset early senescence counter"
        self.tEarlySen[cond712] = 0

        # "Adjust CCx for effects of withered canopy"
        cond713 = (cond71 & (self.CC > self.CCxW))
        self.CCxW[cond713] = self.CC[cond713]

        """
        Calculate canopy size adjusted for micro-advective effects
        """
        # "Check to ensure potential CC is not slightly lower than actual
        cond8 = self.CC_NS < self.CC
        self.CC_NS[cond8] = self.CC[cond8]

        cond81 = (cond8 & (tCC < CanopyDevEnd))
        self.CCxAct_NS[cond81] = self.CC_NS[cond81]

        # Actual (with water stress)
        self.CCadj[self.GrowingSeason] = ((1.72 * self.CC) - (self.CC ** 2) + (0.3 * (self.CC ** 3)))[self.GrowingSeason]

        # Potential (without water stress)
        self.CCadj_NS[self.GrowingSeason] = ((1.72 * self.CC_NS) - (self.CC_NS ** 2) + (0.3 * (self.CC_NS ** 3)))[self.GrowingSeason]

        # "No canopy outside growing season - set various values to zero"
        self.CC[np.logical_not(self.GrowingSeason)] = 0
        self.CCadj[np.logical_not(self.GrowingSeason)] = 0
        self.CC_NS[np.logical_not(self.GrowingSeason)] = 0
        self.CCadj_NS[np.logical_not(self.GrowingSeason)] = 0
        self.CCxW[np.logical_not(self.GrowingSeason)] = 0
        self.CCxAct[np.logical_not(self.GrowingSeason)] = 0
        self.CCxW_NS[np.logical_not(self.GrowingSeason)] = 0
        self.CCxAct_NS[np.logical_not(self.GrowingSeason)] = 0

    def AOS_EvapLayerWaterContent(self, Wevap, MASK = None):
        """
        Function to get water contents in the evaporation layer
        """

        # TODO: remove crop dimension from calculations

        # If no MASK if provided, include all cells
        if MASK is None:
            MASK = np.fill((self.nRotation, self.nLat, self.nLon), True)
        
        for key, value in Wevap.iteritems():
            Wevap[key][MASK] = 0

        # Add dimensions to dz,dzsum
        dz = (self.soil.dz[:,None,None,None] * np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:])
        dzsum = (self.soil.dzsum[:,None,None,None] * np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:])

        # Add compartment dimension to EvapZ
        EvapZ_comp = self.EvapZ[None,:,:,:] * np.ones((self.soil.nComp))[:,None,None,None]

        # expand soil properties to compartments
        th_s   = self.soil.th_s[self.soil.layerIndex,:]
        th_fc  = self.soil.th_fc[self.soil.layerIndex,:]
        th_wp  = self.soil.th_wp[self.soil.layerIndex,:]
        th_dry = self.soil.th_dry[self.soil.layerIndex,:]
        th = self.th

        factor = np.ones((self.soil.nComp, self.nRotation, self.nLat, self.nLon))
        cond1 = (dzsum > EvapZ_comp)
        factor[cond1] = (1 - ((dzsum - EvapZ_comp) / dz))[cond1]
        factor[factor > 1] = 0  # soil layers entirely below EvapZ will have a corresponding factor > 1

        # "Actual water storage in evaporation layer (mm)"
        Wevap['Act'][MASK] = np.sum((factor * 1000 * th * dz), axis=0)[MASK]
        # "Water storage in evaporation layer at saturation (mm)"
        Wevap['Sat'][MASK] = np.sum((factor * 1000 * th_s * dz), axis=0)[MASK]
        # "Water storage in evaporation layer at field capacity (mm)"
        Wevap['Fc'][MASK] = np.sum((factor * 1000 * th_fc * dz), axis=0)[MASK]
        # "Water storage in evaporation layer at wilting point (mm)"
        Wevap['Wp'][MASK] = np.sum((factor * 1000 * th_wp * dz), axis=0)[MASK]
        # "Water storage in evaporation layer at air dry (mm)"
        Wevap['Dry'][MASK] = np.sum((factor * 1000 * th_dry * dz), axis=0)[MASK]

        # Set negative values to zero
        Wevap['Act'] = np.clip(Wevap['Act'], 0, None)

        return Wevap
        
    def AOS_SoilEvaporation(self, meteo, currTimeStep):
        # "Function to calculate daily soil evaporation in AOS"

        # TODO: remove crop dimension from calculations - not required

        # expand soil properties to compartments
        th_dry = self.soil.th_dry[self.soil.layerIndex,:]

        # Add dimensions to dz,dzsum
        dz = (self.soil.dz[:,None,None,None] * np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:])
        dzsum = (self.soil.dzsum[:,None,None,None] * np.ones((self.nRotation, self.nLat, self.nLon))[None,:,:,:])

        arr_zeros = np.zeros((self.nRotation, self.nLat, self.nLon))
        Wevap = {'Act': arr_zeros, 'Sat': arr_zeros, 'Fc': arr_zeros, 'Wp': arr_zeros, 'Dry': arr_zeros}

        # Add crop dimension to some soil parameters...
        # EvapZmin = self.soil.EvapZmin[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None]
        # EvapZmax = self.soil.EvapZmax[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None]
        # REW = self.soil.REW[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None]
        fwcc = self.soil.fwcc[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None]
        Kex = self.soil.Kex[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None]
        # fWrelExp = self.soil.fWrelExp[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None]
        # fevap  = self.soil.fevap[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None]
        
        # ...and some variables
        # Infl = self.Infl[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None]
        SurfaceStorage = self.SurfaceStorage[None,:,:,:] * np.ones((self.crop.nCrop))[:,None,None,None]

        # Add rotation dimension to some crop parameters
        Senescence = self.crop.Senescence[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        
        # Add rotation dimension to meteo vars
        et0 = meteo.referencePotET[None,:,:] * np.ones((self.nRotation))[:,None,None]
        prec = meteo.precipitation[None,:,:] * np.ones((self.nRotation))[:,None,None]
        
        """
        "Prepare stage 2 evaporation (REW gone)"
        """
        # "Only do this if it is first day of the simulation, or if it is first
        # day of growing season and not simulating off-season"
        cond1 = (currTimeStep.timeStepPCR == 1) | ((np.any(self.DAP == 1, axis=0)) & (self.OffSeason == False))
        # "Reset storage in surface soil layer to zero"
        self.Wsurf[cond1] = 0
        # "Set evaporation depth to minimum"
        self.EvapZ[cond1] = self.soil.EvapZmin[cond1]
        # "Trigger stage 2 evaporation"
        self.Stage2[cond1] = True
        # "Get relative water content for start of stage 2 evaporation"
        Wevap = self.AOS_EvapLayerWaterContent(Wevap, MASK = cond1)
        self.Wstage2[cond1] = ((Wevap['Act'] - (Wevap['Fc'] - self.soil.REW)) / (Wevap['Sat'] - (Wevap['Fc'] - self.soil.REW)))[cond1]
        self.Wstage2[cond1] = (np.round((self.Wstage2 * 100)) / 100)[cond1]
        self.Wstage2[cond1] = np.clip(self.Wstage2, 0, None)[cond1]


        """
        "Prepare soil evaporation stage 1"
        """
        # "Adjust water in surface evaporation layer for any infiltration"

        # "Only prepare stage one when rainfall occurs, or when irrigation is
        # triggered (not in net irrigation mode)"

        cond2 = ((prec > 0) | np.any(((self.Irr > 0) & (self.IrrMethod != 4)), axis=0))

        # "Update storage in surface evaporation layer for incoming infiltration"
        cond21 = (cond2 & (self.Infl > 0))
        self.Wsurf[cond21] = self.Infl[cond21]

        # "Water stored in surface evaporation layer cannot exceed REW"
        self.Wsurf[cond21] = np.clip(self.Wsurf, None, self.soil.REW)[cond21]

        # "Reset variables"
        self.Wstage2[cond21] = 0
        self.EvapZ[cond21] = self.soil.EvapZmin[cond21]
        self.Stage2[cond21] = False

        """
        "Calculate potential soil evaporation rate (mm/day)"
        """

        # NB calculation of potential evaporation rate uses crop-specific
        # parameters
        
        # "Adjust time for any delayed development
        tAdj = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        if self.crop.CalendarType == 1:
            tAdj[self.GrowingSeason] = (self.DAP - self.DelayedCDs)[self.GrowingSeason]
        elif self.crop.CalendarType == 2:
            tAdj[self.GrowingSeason] = (self.GDDcum - self.DelayedGDDs)[self.GrowingSeason]

        # "Calculate maximum potential soil evaporation"
        EsPotMax = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        EsPot = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))        
        EsPotMax[self.GrowingSeason] = (Kex * et0 * (1 - self.CCxW * (fwcc / 100)))[self.GrowingSeason]

        # "Calculate potential soil evaporation (given current canopy cover size)
        EsPot[self.GrowingSeason] = (Kex * (1 - self.CCadj) * et0)[self.GrowingSeason]

        # "Adjust potential soil evaporation for effects of withered canopy"
        cond3 = (self.GrowingSeason & (tAdj > Senescence) & (self.CCxAct > 0))
        mult = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))  # preallocate
        cond31 = (cond3 & (self.CC > (self.CCxAct / 2)))
        cond311 = (cond31 & (self.CC <= self.CCxAct))
        mult[cond311] = ((self.CCxAct - self.CC) / (self.CCxAct / 2))[cond311]
        cond32 = (cond3 & np.logical_not(cond31))
        mult[cond32] = 1

        EsPot[cond3] = (EsPot * (1 - self.CCxAct * (fwcc / 100) * mult))[cond3]
        CCxActAdj = ((1.72 * self.CCxAct) + (self.CCxAct ** 2) - 0.3 * (self.CCxAct ** 3))
        EsPotMin = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        EsPotMin[cond3] = (Kex * (1 - CCxActAdj) * et0)[cond3]
        EsPotMin[cond3] = np.clip(EsPotMin, 0, None)[cond3]
        EsPot[cond3] = np.clip(EsPot, EsPotMin, EsPotMax)[cond3]

        cond4 = (self.GrowingSeason & self.PrematSenes)
        EsPot[cond4] = np.clip(EsPot, None, EsPotMax)[cond4]

        # "No canopy cover outside of growing season so potential soil evaporation
        # only depends on reference evapotranspiration"
        EsPot[np.logical_not(self.GrowingSeason)] = (Kex * et0)[np.logical_not(self.GrowingSeason)]

        """
        "Adjust potential soil evaporation for mulches and/or partial wetting"
        """
        # "Mulches"
        EsPotMul = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        cond5 = (SurfaceStorage < 0.000001)

        # "No mulches present"
        cond51 = (cond5 & (self.Mulches == 0))
        EsPotMul[cond51] = EsPot[cond51]

        # "Mulches present (percentage soil surface covered may vary depending
        # on whether within or outside growing season"
        cond52 = (cond5 & (self.Mulches == 1))
        cond521 = (cond52 & self.GrowingSeason)
        EsPotMul[cond521] = (EsPot * (1 - self.fMulch * (self.MulchPctGS / 100)))[cond521]
        cond522 = (cond52 & np.logical_not(self.GrowingSeason))
        EsPotMul[cond522] = (EsPot * (1 - self.fMulch * (self.MulchPctOS / 100)))[cond521]

        # "Surface is flooded - no adjustment of potential soil evaporation for
        # mulches"
        cond6 = np.logical_not(cond5)
        EsPotMul[cond6] = EsPot[cond6]

        """
        "Partial surface wetting by irrigation"
        """
        EsPotIrr = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        # "Only apply adjustment if irrigation occurs and not in net irrigation
        # mode"
        cond7 = ((self.Irr > 0) & (self.IrrMethod != 4))
        # "No adjustment for partial wetting - assume surface is fully wet"
        cond71 = (cond7 & ((prec > 0) | (SurfaceStorage > 0)))
        EsPotIrr[cond71] = EsPot[cond71]
        cond72 = (cond7 & np.logical_not(cond71))
        EsPotIrr[cond72] = (EsPot * (self.WetSurf / 100))[cond72]
        cond8 = np.logical_not(cond7)
        EsPotIrr[cond8] = EsPot[cond8]

        # "Assign minimum value (mulches and partial wetting don't combine"
        EsPot = np.minimum(EsPotIrr,EsPotMul)

        # Sum EsPot along crop dimension
        EsPot = np.sum(EsPot, axis=0)  # crop,rotation,lat,lon -> rotation,lat,lon

        """
        "Surface evaporation"
        """
        # "Initialise actual evaporation counter"
        EsAct = np.zeros((self.nRotation, self.nLat, self.nLon))
        # "Evaporate surface storage"
        cond9 = (self.SurfaceStorage > 0)
        cond91 = (cond9 & (self.SurfaceStorage > EsPot))
        EsAct[cond91] = EsPot[cond91]
        self.SurfaceStorage[cond91] = (self.SurfaceStorage - EsAct)[cond91]

        # Otherwise, "Surface storage is not sufficient to meet all potential soil
        # evaporation"
        cond92 = (cond9 & np.logical_not(cond91))
        EsAct[cond92] = self.SurfaceStorage[cond92]
        # "Update surface storage, evaporation layer depth, stage"
        self.SurfaceStorage[cond92] = 0
        self.Wsurf[cond92] = self.soil.REW[cond92]
        self.Wstage2[cond92] = 0 # not sure if this is correct
        self.EvapZ[cond92] = self.soil.EvapZmin[cond92]
        self.Stage2[cond92] = False
        
        """
        "Stage 1 evaporation"
        """
        # "Determine total water to be extracted"
        ToExtract = EsPot - EsAct

        # "Determine total water to be extracted in stage one (limited by
        # surface layer water storage)"
        ExtractPotStg1 = np.minimum(ToExtract,self.Wsurf)

        # "Extract water"
        cond10 = (ExtractPotStg1 > 0)
        
        # "Find soil compartments covered by evaporation layer"
        comp_sto = np.sum(dzsum < self.soil.EvapZmin, axis=0) + 1
        comp = 0
        while np.any((comp < comp_sto) & (ExtractPotStg1 > 0)):
            
            cond101 = ((comp < comp_sto) & (ExtractPotStg1 > 0))
            
            # "Determine proportion of compartment in evaporation layer"
            factor = np.ones((self.nRotation, self.nLat, self.nLon))
            cond1011 = (cond101 & (dzsum[comp,:] > self.soil.EvapZmin))
            factor[cond1011] = (1 - ((dzsum[comp,:] - self.soil.EvapZmin) / dz[comp,:]))[cond1011]
            cond1012 = (cond101 & np.logical_not(cond1011))
            factor[comp > comp_sto] = 0

            # Water storage at air dry (Wdry) and available water (W) (mm)
            Wdry = np.zeros((self.nRotation, self.nLat, self.nLon))
            W = np.zeros((self.nRotation, self.nLat, self.nLon))
            Wdry = 1000 * th_dry[comp,:] * dz[comp,:]  
            W = 1000 * self.th[comp,:] * dz[comp,:]

            # "Water available in compartment for extraction (mm)"
            AvW = np.zeros((self.nRotation, self.nLat, self.nLon))
            AvW[cond101] = ((W - Wdry) * factor)[cond101]
            AvW = np.clip(AvW, 0, None)

            cond1011 = (cond101 & (AvW >= ExtractPotStg1))
            cond1012 = (cond101 & np.logical_not(cond1011))
            chng = np.zeros((self.nRotation, self.nLat, self.nLon))
            chng[cond1011] = ExtractPotStg1[cond1011]
            chng[cond1012] = AvW[cond1012]

            # Update variables: actual evaporation (EsAct), depth of water in
            # current compartment (W), total water to be extracted (ToExtract)
            # and water to be extracted from surface layer (ExtractPotStg1).
            EsAct += chng       # No need for indices because chng = 0 if conditions 101|102 are not met
            W -= chng
            ToExtract -= chng
            ExtractPotStg1 -= chng

            # "Update water content"
            self.th[comp,:][cond101] = (W / (1000 * dz[comp,:]))[cond101]
            comp += 1           # TODO: check logic

        # "Update surface evaporation layer water balance"
        self.Wsurf[cond10] -= EsAct[cond10]
        cond102 = (cond10 & ((self.Wsurf < 0) | (ExtractPotStg1 > 0.0001)))
        self.Wsurf[cond102] = 0

        # "If surface storage completely depleted, prepare stage 2"
        cond103 = (cond10 & (self.Wsurf < 0.0001))
        # "Get water contents"
        Wevap = self.AOS_EvapLayerWaterContent(Wevap, MASK = cond103)
        self.Wstage2[cond103] = ((Wevap['Act'] - (Wevap['Fc'] - self.soil.REW)) / (Wevap['Sat'] - (Wevap['Fc'] - self.soil.REW)))[cond103]
        self.Wstage2[cond103] = (np.round((self.Wstage2 * 100)) / 100)[cond103]
        self.Wstage2[cond103] = np.clip(self.Wstage2, 0, None)[cond103]

        """
        "Stage 2 evaporation"
        """
        # "Extract water"
        cond11 = (ToExtract > 0)

        if np.any(cond11):

            # "Start stage 2"
            self.Stage2[cond11] = True

            # "Get sub-daily evaporative demand"
            self.EvapTimeSteps = 20  # TODO: add to options???
            Edt = ToExtract / self.EvapTimeSteps

            # "Loop sub-daily time steps"
            for jj in range(self.EvapTimeSteps):
                # "Get current water storage"
                Wevap = self.AOS_EvapLayerWaterContent(Wevap, MASK = cond11)
                # "Get water storage (mm) at start of stage 2 evaporation"
                Wupper = (self.Wstage2 * (Wevap['Sat'] - (Wevap['Fc'] - self.soil.REW)) + (Wevap['Fc'] - self.soil.REW))
                # "Get water storage (mm) when there is no evaporation"
                Wlower = Wevap['Dry']
                # "Get relative depletion of evaporation storage in stage 2"
                Wrel = ((Wevap['Act'] - Wlower) / (Wupper - Wlower))
                # "Check if need to expand evaporative layer"
                cond111 = (cond11 & (self.soil.EvapZmax > self.soil.EvapZmin))
                Wcheck = self.soil.fWrelExp * ((self.soil.EvapZmax - self.EvapZ) / (self.soil.EvapZmax - self.soil.EvapZmin))
                while np.any(cond111 & (Wrel < Wcheck) & (self.EvapZ < self.soil.EvapZmax)):
                    cond1111 = (cond111 & (Wrel < Wcheck) & (self.EvapZ < self.soil.EvapZmax))
                    self.EvapZ[cond111] += 0.001  # "Expand evaporation layer by 1mm"
                    Wevap = AOS_EvapLayerWaterContent(Wevap)  # No need for MASK - only change is EvapZ += 0.001
                    Wupper = (self.Wstage2 * (Wevap['Sat'] - (Wevap['Fc'] - self.soil.REW)) + (Wevap['Fc'] - self.soil.REW))
                    Wlower = Wevap['Dry']
                    Wrel = ((Wevap['Act'] - Wlower) / (Wupper - Wlower))
                    Wcheck = self.soil.fWrelExp * ((self.soil.EvapZmax - self.EvapZ) / (self.soil.EvapZmax - self.soil.EvapZmin))

                # "Get stage 2 evaporation reduction coefficient"
                Kr = ((np.exp(self.soil.fevap * Wrel) - 1) / (np.exp(self.soil.fevap) - 1))
                Kr = np.clip(Kr, None, 1)

                # "Get water to extract"
                ToExtractStg2 = (Kr * Edt)  # Edt is zero in cells which do not need stage 2, so no need for index

                # "Extract water from compartments"
                comp_sto = np.sum(dzsum < self.EvapZ, axis=0) + 1
                comp = 0

                while np.any(cond11 & (comp < comp_sto) & (ToExtractStg2 > 0)):

                    cond111 = (cond11 & (comp < comp_sto) & (ToExtractStg2 > 0))
                    
                    # "Determine proportion of compartment in evaporation layer"
                    factor = np.ones((self.nRotation, self.nLat, self.nLon))
                    cond1111 = (cond111 & (dzsum[comp,:] > self.EvapZ))
                    factor[cond1111] = (1 - ((dzsum[comp,:] - self.EvapZ) / dz[comp,:]))[cond1111]
                    cond1112 = (cond111 & np.logical_not(cond1111))
                    factor[comp > comp_sto] = 0

                    # Water storage at air dry (Wdry) and available water (W) (mm)
                    Wdry = np.zeros((self.nRotation, self.nLat, self.nLon))
                    W = np.zeros((self.nRotation, self.nLat, self.nLon))
                    Wdry = 1000 * th_dry[comp,:] * dz[comp,:]  
                    W = 1000 * self.th[comp,:] * dz[comp,:]

                    # "Water available in compartment for extraction (mm)"
                    AvW = np.zeros((self.nRotation, self.nLat, self.nLon))
                    AvW[cond111] = ((W - Wdry) * factor)[cond111]
                    AvW = np.clip(AvW, 0, None)

                    cond1111 = (cond111 & (AvW >= ExtractPotStg1))
                    cond1112 = (cond111 & np.logical_not(cond1111))
                    chng = np.zeros((self.nRotation, self.nLat, self.nLon))
                    chng[cond1111] = ExtractPotStg1[cond1111]
                    chng[cond1112] = AvW[cond1112]

                    # Update variables: actual evaporation (EsAct), depth of water in
                    # current compartment (W), total water to be extracted (ToExtract)
                    # and water to be extracted from surface layer (ExtractPotStg1).
                    EsAct += chng       # No need for indices because chng = 0 if conditions 101|102 are not met
                    W -= chng
                    ToExtract -= chng
                    ExtractPotStg1 -= chng

                    # "Update water content"
                    self.th[comp,:][cond111] = (W / (1000 * dz[comp,:]))[cond111]
                    comp += 1   # "Update compartment counter"
                
        """
        "Store potential evaporation for irrigation calculations on next day"
        """
        self.Epot = EsPot

    def AOS_AerationStress(self):
        """
        "Function to calculate aeration stress coefficient"
        """

        # Add rotation dimension to LagAer
        LagAer = self.crop.LagAer[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        
        # "Determine aeration stress (root zone)"
        cond1 = (self.thRZ_Act > thRZ_Aer)
        # "Calculate aeration stress coefficient"

        # Preallocate aeration stress coefficient array
        self.Ksa_Aer = np.ones((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))

        cond11 = (cond1 & (self.AerDays < LagAer))
        cond12 = (cond1 & np.logical_not(cond11))
        stress = (1 - ((self.thRZ_Sat - self.thRZ_Act) / (self.thRZ_Sat - self.thRZ_Act)))
        self.Ksa_Aer[cond11] = (1 - ((self.AerDays / 3) * stress))[cond11]
        self.Ksa_Aer[cond12] = ((self.thRZ_Sat - self.thRZ_Act) / (self.thRZ_Sat - self.thRZ_Aer))[cond12]

        # "Increment aeration days counter"
        self.AerDays[cond1] += 1
        cond13 = (cond1 & (self.AerDays > LagAer))
        self.AerDays[cond13] = LagAer[cond13]

        # Otherwise there is no stress, so reset aeration days counter
        self.AerDays[np.logical_not(cond1)] = 0

    # def get_crop_index(self):
    #     # Function to get parameter value of crop for specific growing season

    #     # Dimensions of growing_season are crop,rotation,lat,lon
    #     m,n,p = self.GrowingSeason.shape[1:]
    #     I,J,K = np.ogrid[:m,:n,:p]

    #     crop_index = np.arange(0, self.GrowingSeason.shape[0])[:,None,None,None] * np.ones((m,n,p))[None,:,:,:]
    #     crop_index *= self.GrowingSeason
    #     crop_index = np.max(crop_index, axis=0)
    #     growing_season_index = np.any(self.GrowingSeason, axis=0)
    #     return growing_season_index,crop_index    
        
    def AOS_Transpiration(self, meteo):
        """
        "Function to calculate crop transpiration on current day"
        """

        # Add rotation dimension to some crop parameters, then retrieve value for
        # each rotation (for the moment, we assume that only crop can be grown at
        # a time (i.e. no intercropping). In the future we could account for
        # intercropping by taking a weighted average of the parameters.)

        # Dimensions of growing_season are crop,rotation,lat,lon
        m,n,p = self.GrowingSeason.shape[1:]
        I,J,K = np.ogrid[:m,:n,:p]

        crop_index = np.arange(0, self.GrowingSeason.shape[0])[:,None,None,None] * np.ones((m,n,p))[None,:,:,:]
        crop_index *= self.GrowingSeason
        crop_index = np.max(crop_index, axis=0)
        crop_index = crop_index.astype(int)
        growing_season_index = np.any(self.GrowingSeason, axis=0)
        
        Kcb = (self.crop.Kcb[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[crop_index,I,J,K]
        fage = (self.crop.fage[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[crop_index,I,J,K]
        MaxCanopyCD = (self.crop.MaxCanopyCD[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[crop_index,I,J,K]
        a_Tr = (self.crop.a_Tr[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[crop_index,I,J,K]
        LagAer = (self.crop.LagAer[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None])[crop_index,I,J,K]

        # ***NB values of above parameters for rotations where no crops are
        #    growing equal the value of the first crop***
        
        # Add rotation dimension to ET0
        et0 = meteo.referencePotET[None,:,:] * np.ones((self.nRotation))[:,None,None]

        # """
        # "Calculate transpiration (if in growing season)"
        # """

        # # "1. No prior water stress"

        # # "Update ageing days counter"
        # cond1 = (self.GrowingSeason & (self.DAP > MaxCanopyCD))
        # self.AgeDays_NS = self.DAP - MaxCanopyCD

        # # "Update crop coefficient for ageing of canopy"
        # cond2 = (self.GrowingSeason & (self.AgeDays_NS > 5))
        # cond3 = (self.GrowingSeason & np.logical_not(cond2))
        # Kcb_NS = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        # Kcb_NS[cond2] = (Kcb - ((self.AgeDays_NS - 5) * (fage / 100)) * self.CCxW_NS)[cond2]
        # Kcb_NS[cond3] = Kcb[cond3]

        # # "Update crop coefficient for C02 concentration" ***TODO***
        # # cond4 = (self.GrowingSeason & (CurrentConc > RefConc))
        # # Kcb_NS[cond4] *= (1 - 0.05 * (CurrentConc - RefConc) / (550 - RefConc))[cond4]

        # # "Determine potential transpiration rate (no water stress)"
        # TrPot_NS = Kcb_NS * (self.CCadj_NS) * et0 * self.GrowingSeason

        # # "2. Potential prior water stress and/or delayed development"

        # # "Update ageing days counter"
        # DAPadj = self.DAP - self.DelayedCDs
        # cond5 = (self.GrowingSeason & (DAPadj > MaxCanopyCD))
        # self.AgeDays[cond5] = (DAPadj - MaxCanopyCD)[cond5]

        # # "Update crop coefficient for ageing of canopy"
        # cond6 = (self.GrowingSeason & (self.AgeDays > 5))
        # cond7 = (self.GrowingSeason & np.logical_not(cond6))
        # Kcb[cond6] = (Kcb - ((self.AgeDays - 5) * (fage / 100)) * self.CCxW)[cond6]
        # Kcb[cond7] = Kcb[cond7]

        # # "Update crop coefficient for CO2 concentration" ***TODO***
        # # cond8 = (self.GrowingSeason & (CurrentConc > RefConc))
        # # Kcb[cond8] *= (1 - 0.05 * (CurrentConc - RefConc) / (550 - RefConc))[cond4][cond8]

        # # "Determine potential transpiration rate"
        # TrPot0 = Kcb * (self.CCadj) * et0 * self.GrowingSeason

        # # "Correct potential transpiration for dying green canopy effects"
        # cond9 = (self.GrowingSeason & (self.CC < self.CCxW))
        # cond91 = (cond9 & (self.CCxW > 0.001) & (self.CC > 0.001))
        # TrPot0[cond91] *= ((self.CC / self.CCxW) ** a_Tr)[cond91]

        # """
        # "Calculate surface layer transpiration"
        # """
        # cond10 = (self.GrowingSeason & (SurfaceStorage > 0) & (self.DaySubmerged < LagAer))
        # # Initialise variables
        # TrAct0 = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        # TrPot = np.zeros((self.crop.nCrop, self.nRotation, self.nLat, self.nLon))
        
        # # "Update submergence days counter"
        # print cond10.shape
        # print self.DaySubmerged.shape
        # self.DaySubmerged[cond10] += 1

        # # "Update anaerobic conditions counter for each compartment"
        # cond10_comp = np.broadcast_to(cond10, self.AerDaysComp.shape)  # Add compartment dimension
        # self.AerDaysComp[cond10_comp] += 1 
        # LagAer_comp = np.broadcast_to(LagAer, self.AerDaysComp.shape)
        # self.AerDaysComp[cond10_comp] = np.clip(self.AerDaysComp, None, LagAer_comp)[cond10_comp]

        # # "Reduce actual transpiration that is possible to account for aeration
        # # stress due to extended submergence"
        # fSub = 1 - (self.DaySubmerged / LagAer)

        # # "Transpiration occurs from surface storage
        # cond101 = (cond10 & (SurfaceStorage > (fSub * TrPot0)))
        # SurfaceStorage[cond101] -= (fSub * TrPot0)[cond101]
        # self.SurfaceStorage = np.max(SurfaceStorage, axis=0)  # TODO: check this works as expected
        # TrAct0[cond101] = (fSub * TrPot0)[cond101]

        # # Otherwise, "No transpiration from surface storage"
        # cond102 = (cond10 & np.logical_not(cond101))
        # TrAct0[cond102] = 0

        # # "More water can be extracted from soil profile for transpiration"
        # cond103 = (cond10 & (TrAct0 < (fSub * TrPot0)))
        # TrPot[cond103] = ((fSub * TrPot0) - TrAct0)[cond103]

        # """
        # "Update potential root zone transpiration for water stress"
        # """
        # # "Determine root zone water content"
        # self.RootZoneWater()    # TODO: should we change RootZoneWater to return variables rather than update object?

        # # "Calculate water stress coefficients"
        # self.AOS_WaterStress(meteo, beta = True)

        # # "Calculate aeration stress coefficients"
        # self.AOS_AerationStress()

        # # "Maximum stress effect"
        # Ks = np.minimum(self.Ksw_StoLin, self.Ksa_Aer)

        # # "Update potential transpiration in root zone"
        # cond11 = (self.GrowingSeason & (self.IrrMethod != 4))
        # TrPot[cond11] *= Ks[cond11]

        # """
        # "Determine compartments covered by root zone"
        # """
        # # TODO
        
        
    def update(self, meteo, groundwater, currTimeStep):

        """
        Update season counter
        """
        
        # TODO: consider whether CropSequence should be cast to nLat, nLon upon loading

        # Broadcast to include rotation dimension, then multiply by crop sequence
        # to ensure the counter is only updated in grid cells where the crop is
        # grown (TODO: check this is performing as expected)
        cond1 = ((currTimeStep.doy == self.crop.PlantingDate)[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]) *\
                (self.CropSequence[:,:,None,None] * np.ones((self.nLat,self.nLon))[None,None,:,:])
        cond1 = cond1.astype(bool)
        self.SeasonCounter[cond1] += 1

        """
        Update growing season
        """

        cond1 = (self.SeasonCounter > 0)

        # Broadcast to include rotation dimension, coerce to boolean 
        cond2 = ((self.crop.PlantingDate <= currTimeStep.doy) & (currTimeStep.doy <= self.crop.HarvestDate))
        cond2 = cond2[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        cond2 = cond2.astype(bool)

        cond3 = (self.crop.PlantingDate > self.crop.HarvestDate)
        cond3 = cond3[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        cond3 = cond3.astype(bool)

        cond4 = ((self.crop.PlantingDate <= currTimeStep.doy) | (currTimeStep.doy <= self.crop.HarvestDate))
        cond4 = cond4[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]
        cond4 = cond4.astype(bool)

        self.GrowingSeason = ((cond1 & cond2) | (cond1 & cond3 & cond4))  # crop,rotation,lat,lon 

        """
        Increment time counters
        """
        
        # update DAP (day after planting)
        self.DAP[self.GrowingSeason] += 1
        self.DAP[np.logical_not(self.GrowingSeason)] = 0

        # Compute growing degree day
        self.crop.GrowingDegreeDay(meteo, currTimeStep)  # crop,nlat,nlon
        gdd = self.crop.GDD[:,None,:,:] * np.ones((self.nRotation))[None,:,None,None]

        self.GDDcum[self.GrowingSeason] += gdd[self.GrowingSeason]
        self.GDDcum[np.logical_not(self.GrowingSeason)] = 0

        # TODO: only do this if a given option is set
        self.AOS_CheckGroundwaterTable(groundwater, currTimeStep)
        self.PreIrrigation()
        self.Drainage()
        self.RainfallPartition(meteo)
        self.Irrigation(meteo)
        self.Infiltration()
        self.CapillaryRise(groundwater)
        self.AOS_Germination()
        self.AOS_GrowthStage()
        self.AOS_RootDevelopment(groundwater)
        self.AOS_CanopyCover(meteo)
        self.AOS_SoilEvaporation(meteo, currTimeStep)
        self.AOS_Transpiration(meteo)
