#!/usr/bin/env python
# -*- coding: utf-8 -*-

netcdf_short_name = {}
netcdf_dimensions = {}
netcdf_unit       = {}
# netcdf_monthly_total_unit = {} 
# netcdf_yearly_total_unit  = {}
netcdf_long_name  = {}
description       = {}
comment           = {}
latex_symbol      = {}

aquacrop_variable_name = 'precipitation'
netcdf_short_name[aquacrop_variable_name] = 'precipitation'
netcdf_dimensions[aquacrop_variable_name] = ('time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'm'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'th'
netcdf_short_name[aquacrop_variable_name] = 'water_content'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','depth','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Wr'
netcdf_short_name[aquacrop_variable_name] = 'root_zone_water_content'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'zGW'
netcdf_short_name[aquacrop_variable_name] = 'groundwater_depth'
netcdf_dimensions[aquacrop_variable_name] = ('time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'm'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'SurfaceStorage'
netcdf_short_name[aquacrop_variable_name] = 'surface_layer_water_content'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Irr'
netcdf_short_name[aquacrop_variable_name] = 'irrigation_depth'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Infl'
netcdf_short_name[aquacrop_variable_name] = 'infiltration'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Runoff'
netcdf_short_name[aquacrop_variable_name] = 'surface_runoff'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'DeepPerc'
netcdf_short_name[aquacrop_variable_name] = 'deep_percolation'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'CrTot'
netcdf_short_name[aquacrop_variable_name] = 'capillary_rise'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'GwIn'
netcdf_short_name[aquacrop_variable_name] = 'groundwater_inflow'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'EsAct'
netcdf_short_name[aquacrop_variable_name] = 'evaporation'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Epot'
netcdf_short_name[aquacrop_variable_name] = 'potential_evaporation'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'TrAct'
netcdf_short_name[aquacrop_variable_name] = 'transpiration'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Tpot'
netcdf_short_name[aquacrop_variable_name] = 'potential_transpiration'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

# crop growth
# ###########

aquacrop_variable_name = 'GDD'
netcdf_short_name[aquacrop_variable_name] = 'growing_degree_days'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'days'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'GDDcum'
netcdf_short_name[aquacrop_variable_name] = 'cumulative_growing_degree_days'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'days'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Zroot'
netcdf_short_name[aquacrop_variable_name] = 'root_depth'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'm'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'CC'
netcdf_short_name[aquacrop_variable_name] = 'canopy_cover'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'CC_NS'
netcdf_short_name[aquacrop_variable_name] = 'reference_canopy_cover'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'B'
netcdf_short_name[aquacrop_variable_name] = 'biomass'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 kg m-2'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'B_NS'
netcdf_short_name[aquacrop_variable_name] = 'reference_biomass'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 kg m-2'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'HI'
netcdf_short_name[aquacrop_variable_name] = 'harvest_index'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'HIadj'
netcdf_short_name[aquacrop_variable_name] = 'adjusted_harvest_index'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Y'
netcdf_short_name[aquacrop_variable_name] = 'crop_yield'
netcdf_dimensions[aquacrop_variable_name] = ('crop','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 kg 1e4 m-2'  # tonne/hectare
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

# not used:

# # dump variables:

# # Counters
# # ########

# aquacrop_variable_name = 'AgeDays'
# netcdf_short_name[aquacrop_variable_name] = 'canopy_age'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'AgeDays_NS'
# netcdf_short_name[aquacrop_variable_name] = 'canopy_age_no_stress'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # AerDays
# aquacrop_variable_name = 'AerDays'
# netcdf_short_name[aquacrop_variable_name] = 'aeration_counter'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # IrrCum
# aquacrop_variable_name = 'IrrCum'
# netcdf_short_name[aquacrop_variable_name] = 'cumulative_irrigation'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # IrrNetCum
# aquacrop_variable_name = 'IrrNetCum'
# netcdf_short_name[aquacrop_variable_name] = 'cumulative_net_irrigation'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # DelayedGDDs
# aquacrop_variable_name = 'DelayedGDDs'
# netcdf_short_name[aquacrop_variable_name] = 'delayed_growth_growing_degree_days'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'DelayedCDs'
# netcdf_short_name[aquacrop_variable_name] = 'delayed_growth_calendar_days'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'PctLagPhase'
# netcdf_short_name[aquacrop_variable_name] = 'lag_phase_progression'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'percent'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'tEarlySen'
# netcdf_short_name[aquacrop_variable_name] = 'early_senescence_counter'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # GDDcum
# aquacrop_variable_name = 'GDDcum'
# netcdf_short_name[aquacrop_variable_name] = 'cumulative_growing_degree_days'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'DaySubmerged'
# netcdf_short_name[aquacrop_variable_name] = 'time_submerged'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'DAP'
# netcdf_short_name[aquacrop_variable_name] = 'days_after_planting'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Epot
# aquacrop_variable_name = 'Epot'
# netcdf_short_name[aquacrop_variable_name] = 'potential_evaporation'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'mm'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Tpot
# aquacrop_variable_name = 'Tpot'
# netcdf_short_name[aquacrop_variable_name] = 'potential_transpiration'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None
        
# # States
# # ######

# aquacrop_variable_name = 'WTinSoil'
# netcdf_short_name[aquacrop_variable_name] = 'water_table_in_soil_profile_flag'
# netcdf_dimensions                         = ('crop','depth','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'PreAdj'
# netcdf_short_name[aquacrop_variable_name] = 'harvest_index_adjusted_for_water_stress_flag'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # CropMature
# aquacrop_variable_name = 'CropMature'
# netcdf_short_name[aquacrop_variable_name] = 'crop_mature_flag'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'CropDead'
# netcdf_short_name[aquacrop_variable_name] = 'crop_dead_flag'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Germination
# aquacrop_variable_name = 'Germination'
# netcdf_short_name[aquacrop_variable_name] = 'crop_germinated_flag'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # PrematSenes
# aquacrop_variable_name = 'PrematSenes'
# netcdf_short_name[aquacrop_variable_name] = 'premature_senescence_flag'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # HarvestFlag
# aquacrop_variable_name = 'HarvestFlag'
# netcdf_short_name[aquacrop_variable_name] = 'harvest_flag'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None


# # Harvest index
# # *************
# # Stage
# aquacrop_variable_name = 'Stage'
# netcdf_short_name[aquacrop_variable_name] = 'crop_growth_stage'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'Fpre'
# netcdf_short_name[aquacrop_variable_name] = 'pre_anthesis_water_stress_adjustment_factor'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'Fpost'
# netcdf_short_name[aquacrop_variable_name] = 'post_anthesis_water_stress_adjustment_factor'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Equation 3.12m in AquaCrop reference manual
# aquacrop_variable_name = 'fpost_dwn'
# netcdf_short_name[aquacrop_variable_name] = 'post_anthesis_water_stress_adjustment_factor_negative'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Equation 3.12l in AquaCrop reference manual
# aquacrop_variable_name = 'fpost_upp'
# netcdf_short_name[aquacrop_variable_name] = 'post_anthesis_water_stress_adjustment_factor_positive'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # aquacrop_variable_name = 'HIcor_Asum'
# # netcdf_short_name[aquacrop_variable_name] = 'HIcor_Asum'
# # netcdf_dimensions                         = ('crop','lat','lon')
# # netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# # netcdf_long_name[aquacrop_variable_name]  = None
# # description[aquacrop_variable_name]       = None
# # comment[aquacrop_variable_name]           = None
# # latex_symbol[aquacrop_variable_name]      = None

# # aquacrop_variable_name = 'HIcor_Bsum'
# # netcdf_short_name[aquacrop_variable_name] = 'HIcor_Bsum'
# # netcdf_dimensions                         = ('crop','lat','lon')
# # netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# # netcdf_long_name[aquacrop_variable_name]  = None
# # description[aquacrop_variable_name]       = None
# # comment[aquacrop_variable_name]           = None
# # latex_symbol[aquacrop_variable_name]      = None

# # Fpol
# aquacrop_variable_name = 'Fpol'
# netcdf_short_name[aquacrop_variable_name] = 'pollination_failure_adjustment_factor'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'sCor1'
# netcdf_short_name[aquacrop_variable_name] = 'harvest_index_leaf_expansion_adjustment'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'sCor2'
# netcdf_short_name[aquacrop_variable_name] = 'harvest_index_stomatal_closure_adjustment'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'GrowthStage'
# netcdf_short_name[aquacrop_variable_name] = 'crop_growth_stage'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'AerDaysComp'
# netcdf_short_name[aquacrop_variable_name] = 'compartment_aeration_counter'
# netcdf_dimensions                         = ('crop','depth','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'd'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'TrRatio'
# netcdf_short_name[aquacrop_variable_name] = 'transpiration_ratio'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'CC'
# netcdf_short_name[aquacrop_variable_name] = 'canopy_cover'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Actual canopy cover adjusted for micro-advective effects
# aquacrop_variable_name = 'CCadj'
# netcdf_short_name[aquacrop_variable_name] = 'adjusted_canopy_cover'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Reference canopy cover (i.e. canopy cover with no water stress)
# aquacrop_variable_name = 'CC_NS'
# netcdf_short_name[aquacrop_variable_name] = 'reference_canopy_cover'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Reference canopy cover adjusted for micro-advective effects
# aquacrop_variable_name = 'CCadj_NS'
# netcdf_short_name[aquacrop_variable_name] = 'adjusted_reference_canopy_cover'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # B
# aquacrop_variable_name = 'B'
# netcdf_short_name[aquacrop_variable_name] = 'biomass'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1e-3 kg m-2'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # B_NS
# aquacrop_variable_name = 'B_NS'
# netcdf_short_name[aquacrop_variable_name] = 'reference_biomass'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1e-3 kg m-2'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Harvest index
# aquacrop_variable_name = 'HI'
# netcdf_short_name[aquacrop_variable_name] = 'harvest_index'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Harvest index adjusted for stress effects
# aquacrop_variable_name = 'HIadj'
# netcdf_short_name[aquacrop_variable_name] = 'adjusted_harvest_index'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Maximum canopy cover
# aquacrop_variable_name = 'CCxAct'
# netcdf_short_name[aquacrop_variable_name] = 'maximum_canopy_cover'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Reference maximum canopy cover
# aquacrop_variable_name = 'CCxAct_NS'
# netcdf_short_name[aquacrop_variable_name] = 'reference_maximum_canopy_cover'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'CCxW'
# netcdf_short_name[aquacrop_variable_name] = 'withered_canopy_cover'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# aquacrop_variable_name = 'CCxW_NS'
# netcdf_short_name[aquacrop_variable_name] = 'reference_withered_canopy_cover'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # CCxEarlySen
# aquacrop_variable_name = 'CCxEarlySen'
# netcdf_short_name[aquacrop_variable_name] = 'maximum_canopy_cover_early_senescence'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # aquacrop_variable_name = 'CCprev'
# # netcdf_short_name[aquacrop_variable_name] = 'water_content'
# # netcdf_dimensions                         = ('crop','depth','lat','lon')
# # netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# # netcdf_long_name[aquacrop_variable_name]  = None
# # description[aquacrop_variable_name]       = None
# # comment[aquacrop_variable_name]           = None
# # latex_symbol[aquacrop_variable_name]      = None

# # correction factor
# aquacrop_variable_name = 'rCor'
# netcdf_short_name[aquacrop_variable_name] = 'water_extraction_correction_factor'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Zroot
# aquacrop_variable_name = 'Zroot'
# netcdf_short_name[aquacrop_variable_name] = 'effective_rooting_depth'
# netcdf_dimensions                         = ('crop','depth','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # CC0adj
# aquacrop_variable_name = 'CC0adj'
# netcdf_short_name[aquacrop_variable_name] = 'adjusted_canopy_cover_at_emergence'
# netcdf_dimensions                         = ('crop','depth','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # depth of water stored at the surface
# aquacrop_variable_name = 'SurfaceStorage'
# netcdf_short_name[aquacrop_variable_name] = 'surface_storage'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # initial depth of water at the surface
# aquacrop_variable_name = 'SurfaceStorageIni'
# netcdf_short_name[aquacrop_variable_name] = 'initial_surface_storage'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # th
# aquacrop_variable_name = 'th'
# netcdf_short_name[aquacrop_variable_name] = 'water_content'
# netcdf_dimensions                         = ('crop','depth','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None


# # Wsurf
# aquacrop_variable_name = 'Wsurf'
# netcdf_short_name[aquacrop_variable_name] = 'surface_evaporation_layer_storage'
# netcdf_dimensions                         = ('crop','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # EvapZ
# aquacrop_variable_name = 'EvapZ'
# netcdf_short_name[aquacrop_variable_name] = 'evaporation_layer_depth'
# netcdf_dimensions                         = ('crop','depth','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Flag to indicate stage 2 soil evaporation calculation procedure has been initiated
# aquacrop_variable_name = 'Stage2'
# netcdf_short_name[aquacrop_variable_name] = 'soil_evaporation_stage_2_flag'
# netcdf_dimensions                         = ('crop','depth','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = '1'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None

# # Relative water content of evaporation layer at start of stage 2 soil evaporation
# aquacrop_variable_name = 'Wstage2'
# netcdf_short_name[aquacrop_variable_name] = 'soil_evaporation_stage_2_initial_water_content'
# netcdf_dimensions                         = ('crop','depth','lat','lon')
# netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
# netcdf_long_name[aquacrop_variable_name]  = None
# description[aquacrop_variable_name]       = None
# comment[aquacrop_variable_name]           = None
# latex_symbol[aquacrop_variable_name]      = None
