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

aquacrop_variable_name = 'th'
netcdf_short_name[aquacrop_variable_name] = 'water_content'
netcdf_dimensions                         = ('rotation','time','depth','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'm3 m-3'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Wr'
netcdf_short_name[aquacrop_variable_name] = 'root_zone_water_content'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'zGW'
netcdf_short_name[aquacrop_variable_name] = 'groundwater_depth'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'm'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'SurfaceStorage'
netcdf_short_name[aquacrop_variable_name] = 'surface_layer_water_content'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Irr'
netcdf_short_name[aquacrop_variable_name] = 'irrigation_depth'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Infl'
netcdf_short_name[aquacrop_variable_name] = 'infiltration'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Runoff'
netcdf_short_name[aquacrop_variable_name] = 'surface_runoff'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'DeepPerc'
netcdf_short_name[aquacrop_variable_name] = 'deep_percolation'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'CR'
netcdf_short_name[aquacrop_variable_name] = 'capillary_rise'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'GwIn'
netcdf_short_name[aquacrop_variable_name] = 'groundwater_inflow'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'EsAct'
netcdf_short_name[aquacrop_variable_name] = 'evaporation'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'EsPot'
netcdf_short_name[aquacrop_variable_name] = 'potential_evaporation'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'TrAct'
netcdf_short_name[aquacrop_variable_name] = 'transpiration'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'TrPot'
netcdf_short_name[aquacrop_variable_name] = 'potential_transpiration'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 m'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

# crop growth
# ###########

aquacrop_variable_name = 'GDD'
netcdf_short_name[aquacrop_variable_name] = 'growing_degree_days'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'days'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'GDDcum'
netcdf_short_name[aquacrop_variable_name] = 'cumulative_growing_degree_days'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'days'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Zroot'
netcdf_short_name[aquacrop_variable_name] = 'root_depth'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = 'm'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'CC'
netcdf_short_name[aquacrop_variable_name] = 'canopy_cover'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'CC_NS'
netcdf_short_name[aquacrop_variable_name] = 'reference_canopy_cover'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'B'
netcdf_short_name[aquacrop_variable_name] = 'biomass'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 kg m-2'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'B_NS'
netcdf_short_name[aquacrop_variable_name] = 'reference_biomass'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 kg m-2'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'HI'
netcdf_short_name[aquacrop_variable_name] = 'harvest_index'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'HIadj'
netcdf_short_name[aquacrop_variable_name] = 'adjusted_harvest_index'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1'
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None

aquacrop_variable_name = 'Y'
netcdf_short_name[aquacrop_variable_name] = 'crop_yield'
netcdf_dimensions                         = ('rotation','time','lat','lon')
netcdf_unit[aquacrop_variable_name]       = '1e-3 kg 1e4 m-2'  # tonne/hectare
netcdf_long_name[aquacrop_variable_name]  = None
description[aquacrop_variable_name]       = None
comment[aquacrop_variable_name]           = None
latex_symbol[aquacrop_variable_name]      = None
