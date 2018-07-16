#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model
import os
import shutil
import sys
import math
import gc

import pcraster as pcr
import VirtualOS as vos

from Model import Model
from BiomassAccumulation import *
from CanopyCover import *
from CapillaryRise import *
from CarbonDioxide import *
from CheckGroundwaterTable import *
from CropParameters import *
from CropYield import *
from Drainage import *
from FieldMgmtParameters import *
from Germination import *
from Groundwater import *
from GrowthStage import *
from HarvestIndex import *
from HIrefCurrentDay import *
from Infiltration import *
from Inflow import *
from InitialCondition import *
from IrrigationMgmtParameters import *
from Irrigation import *
from Meteo import *
from PreIrrigation import *
from RainfallPartition import *
from RootDevelopment import *
from RootZoneWater import *
from SoilAndTopoParameters import *
from SoilEvaporation import *
from Evapotranspiration import *
from TemperatureStress import *
from Transpiration import *
from WaterStress import *

import logging
logger = logging.getLogger(__name__)

class FAO56(Model):

    def initial(self):
        
        self.meteo_module = Meteo(self)
        self.groundwater_module = Groundwater(self)
        self.carbon_dioxide_module = CarbonDioxide(self)

        self.crop_parameters_module = FAO56CropParameters(self)
        self.field_mgmt_parameters_module = FieldMgmtParameters(self)
        self.irrigation_mgmt_parameters_module = IrrigationMgmtParameters(self)
        self.soil_parameters_module = SoilAndTopoParameters(self)
        
        self.initial_condition_module = InitialCondition(self)
        self.check_groundwater_table_module = CheckGroundwaterTable(self)
        self.pre_irrigation_module = PreIrrigation(self)
        self.drainage_module = Drainage(self)
        self.rainfall_partition_module = RainfallPartition(self)
        self.root_zone_water_module = FAO56RootZoneWater(self)
        self.irrigation_module = Irrigation(self)
        self.infiltration_module = Infiltration(self)
        self.capillary_rise_module = CapillaryRise(self)
        # self.germination_module = Germination(self)
        self.growth_stage_module = FAO56GrowthStage(self)
        self.root_development_module = FAO56RootDevelopment(self)
        # self.water_stress_module = WaterStress(self)
        # self.canopy_cover_module = CanopyCover(self)
        # self.soil_evaporation_module = SoilEvaporation(self)
        # self.transpiration_module = Transpiration(self)
        self.evapotranspiration_module = FAO56Evapotranspiration(self)
        self.inflow_module = Inflow(self)
        # self.HI_ref_current_day_module = HIrefCurrentDay(self)
        # self.biomass_accumulation_module = BiomassAccumulation(self)
        # self.temperature_stress_module = TemperatureStress(self)
        # self.harvest_index_module = HarvestIndex(self)
        self.crop_yield_module = FAO56CropYield(self)

        # initialize modules
        self.meteo_module.initial()
        self.groundwater_module.initial()
        self.carbon_dioxide_module.initial()

        self.crop_parameters_module.initial()
        self.field_mgmt_parameters_module.initial()
        self.irrigation_mgmt_parameters_module.initial()
        self.soil_parameters_module.initial()
        
        self.initial_condition_module.initial()
        self.check_groundwater_table_module.initial()
        self.pre_irrigation_module.initial()
        self.drainage_module.initial()
        self.rainfall_partition_module.initial()
        self.root_zone_water_module.initial()
        self.irrigation_module.initial()
        self.infiltration_module.initial()
        self.capillary_rise_module.initial()
        # self.germination_module.initial()
        self.growth_stage_module.initial()
        self.root_development_module.initial()
        # self.water_stress_module.initial()
        # self.canopy_cover_module.initial()
        # self.soil_evaporation_module.initial()
        # self.transpiration_module.initial()
        self.evapotranspiration_module.initial()
        self.inflow_module.initial()
        # self.HI_ref_current_day_module.initial()
        # self.biomass_accumulation_module.initial()
        # self.temperature_stress_module.initial()
        # self.harvest_index_module.initial()
        self.crop_yield_module.initial()
        
    def dynamic(self):
        """Function to update model state for current time step"""

        logger.info("Reading forcings for time %s", self._modelTime)
        self.meteo_module.dynamic()
        self.groundwater_module.dynamic()
        # self.carbon_dioxide_module.dynamic()

        self.crop_parameters_module.dynamic()
        self.irrigation_mgmt_parameters_module.dynamic()

        self.initial_condition_module.dynamic()
        self.check_groundwater_table_module.dynamic()
        self.pre_irrigation_module.dynamic()
        self.drainage_module.dynamic()
        self.rainfall_partition_module.dynamic()

        # ***at this point we need to run the livelihood model, and estimate
        # the maximum potential irrigation the farmer can afford to apply
        # given the current groudwater level***
                
        self.root_zone_water_module.dynamic()        
        self.irrigation_module.dynamic()
        self.infiltration_module.dynamic()
        self.capillary_rise_module.dynamic()

        self.growth_stage_module.dynamic()
        
        # self.root_development_module.dynamic()
        self.root_zone_water_module.dynamic()
        # self.water_stress_module.dynamic(beta=True)

        # self.transpiration_module.dynamic()  # TODO: replace with evapotranspiration module (possibly combine with water stress)
        self.evapotranspiration_module.dynamic()
        self.inflow_module.dynamic()

        # self.HI_ref_current_day_module.dynamic()
        # self.biomass_accumulation_module.dynamic()

        # self.root_zone_water_module.dynamic()
        # self.water_stress_module.dynamic(beta=True)
        # self.temperature_stress_module.dynamic()
        # self.harvest_index_module.dynamic()
        self.crop_yield_module.dynamic()

        self.root_zone_water_module.dynamic()
        self.pre_irrigation_module.add_pre_irrigation()
