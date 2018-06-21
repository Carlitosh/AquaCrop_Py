#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
import shutil
import sys
import math
import gc
import numpy as np

def evap_layer_water_content(th, th_s, th_fc, th_wp, th_dry, dz, evapz):
    """Function to get water contents in the evaporation layer"""

    dims = th.shape
    arr_ones = np.ones((dims[1], dims[2], dims[3]))
    dz = dz[:,None,None,None] * arr_ones
    dzsum = np.cumsum(dz, axis=0)

    # Find compartments covered by evaporation layer
    comp_sto = (np.round((dzsum - dz) * 1000) < np.round(evapz * 1000))
    factor = 1 - ((dzsum - evapz) / dz)
    factor = np.clip(factor, 0, 1) * comp_sto

    # Water storages in evaporation layer (mm)
    Wevap = {}
    Wevap_Act = np.sum((factor * 1000 * th * dz), axis=0)
    Wevap['Act'] = np.clip(Wevap_Act, 0, None)
    Wevap['Sat'] = np.sum((factor * 1000 * th_s * dz), axis=0)
    Wevap['Fc'] = np.sum((factor * 1000 * th_fc * dz), axis=0)
    Wevap['Wp'] = np.sum((factor * 1000 * th_wp * dz), axis=0)
    Wevap['Dry'] = np.sum((factor * 1000 * th_dry * dz), axis=0)
    return Wevap

def prepare_stage_two_evaporation(th, th_s, th_fc, th_wp, th_dry, dz, REW, EvapZ):
    Wevap = evap_layer_water_content(th, th_s, th_fc, th_wp, th_dry, dz, EvapZ)
    Wstage2 = ((Wevap['Act'] - (Wevap['Fc'] - REW)) / (Wevap['Sat'] - (Wevap['Fc'] - REW)))
    Wstage2 = (np.round((Wstage2 * 100)) / 100)
    Wstage2 = np.clip(Wstage2, 0, None)
    return Wstage2

def surface_evaporation(surface_storage, EsPot):

    dims = surface_storage.shape
    nr, nlat, nlon = dims[0], dims[1], dims[2]
    
    # Initialise actual evaporation counter
    EsActSurf = np.zeros((nr, nlat, nlon))
    cond9 = (surface_storage > 0)
    cond91 = (cond9 & (surface_storage > EsPot))
    EsActSurf[cond91] = EsPot[cond91]
    surface_storage[cond91] = (surface_storage - EsActSurf)[cond91]
    cond92 = (cond9 & np.logical_not(cond91))
    EsActSurf[cond92] = surface_storage[cond92]
    return EsActSurf

def extract_water(th, th_dry, dz, dzsum, EvapZmin, EsAct, ToExtract, ToExtractStg):

    # Extract water
    # cond10 = (ToExtractStg > 0)
    dims = th.shape
    nc, nr, nlat, nlon = dims[0], dims[1], dims[2], dims[3]
    thnew = np.copy(th)
    
    # Determine fraction of compartments covered by evaporation layer
    comp_sto = (np.round((dzsum - dz) * 1000) < np.round(EvapZmin * 1000))
    factor = 1 - ((dzsum - EvapZmin) / dz)
    factor = np.clip(factor, 0, 1) * comp_sto

    comp_sto = np.sum(comp_sto, axis=0)
    comp = 0
    while np.any((comp < comp_sto) & (ToExtractStg > 0) & (ToExtract > 0)):

        cond101 = ((comp < comp_sto) & (ToExtractStg > 0) & (ToExtract > 0))

        # Water available in compartment for extraction (mm)
        Wdry = 1000 * th_dry[comp,:] * dz[comp,:]  
        W = 1000 * thnew[comp,:] * dz[comp,:]
        AvW = np.zeros((nr, nlat, nlon))
        AvW[cond101] = ((W - Wdry) * factor[comp,:])[cond101]
        AvW = np.clip(AvW, 0, None)

        # Determine amount by which to adjust variables
        cond1011 = (cond101 & (AvW >= ToExtractStg))
        EsAct[cond1011] += ToExtractStg[cond1011]
        W[cond1011] -= ToExtractStg[cond1011]
        ToExtract[cond1011] -= ToExtractStg[cond1011]
        ToExtractStg[cond1011] = 0

        cond1012 = (cond101 & np.logical_not(cond1011))
        EsAct[cond1012] += AvW[cond1012]
        W[cond1012] -= AvW[cond1012]
        ToExtract[cond1012] -= AvW[cond1012]
        ToExtractStg[cond1012] -= AvW[cond1012]

        # Update water content
        thnew[comp,:][cond101] = (W / (1000 * dz[comp,:]))[cond101]
        comp += 1
        
    return thnew, EsAct, ToExtract, ToExtractStg

def soil_evaporation_stage_one(th, th_s, th_fc, th_wp, th_dry, EsPot, EsAct, dz, dzsum, EvapZmin, REW, EvapZ, Wsurf, Wstage2):

    dims = th.shape
    nc, nr, nlat, nlon = dims[0], dims[1], dims[2], dims[3]
    thnew = np.copy(th)
    dz = dz[:,None,None,None] * np.ones((nr, nlat, nlon))[None,:,:,:]
    dzsum = dzsum[:,None,None,None] * np.ones((nr, nlat, nlon))[None,:,:,:]
        
    # Determine total water to be extracted
    ToExtract = EsPot - EsAct

    # Determine total water to be extracted in stage one (limited by
    # surface layer water storage)
    ExtractPotStg1 = np.minimum(ToExtract, Wsurf)

    # Extract water
    cond10 = (ExtractPotStg1 > 0)

    thnew, EsAct, ToExtract, ExtractPotStg1 = extract_water(
        thnew, th_dry, dz, dzsum, EvapZmin, EsAct, ToExtract, ExtractPotStg1)
    
    # Update surface evaporation layer water balance
    Wsurf[cond10] -= EsAct[cond10]
    cond102 = (cond10 & ((Wsurf < 0) | (ExtractPotStg1 > 0.0001)))
    Wsurf[cond102] = 0

    # If surface storage completely depleted, prepare stage 2 evaporation
    cond103 = (cond10 & (Wsurf < 0.0001))
    Wstage2_tmp = prepare_stage_two_evaporation(thnew, th_s, th_fc, th_wp, th_dry, dz[:,0,0,0], REW, EvapZ)
    Wstage2[cond103] = Wstage2_tmp[cond103]  # TODO: add MASK argument to function?
    return thnew, Wsurf, ToExtract, Wstage2

def relative_depletion(th, th_s, th_fc, th_wp, th_dry, dz, EvapZ, Wstage2, REW):

    # Get current water storage
    Wevap = evap_layer_water_content(th, th_s, th_fc, th_wp, th_dry, dz, EvapZ)
    # Get water storage (mm) at start of stage 2 evaporation
    Wupper = (Wstage2 * (Wevap['Sat'] - (Wevap['Fc'] - REW)) + (Wevap['Fc'] - REW))
    # Get water storage (mm) when there is no evaporation
    Wlower = np.copy(Wevap['Dry'])
    # Get relative depletion of evaporation storage in stage 2
    Wrel_divd = (Wevap['Act'] - Wlower)
    Wrel_divs = (Wupper - Wlower)
    Wrel = np.divide(Wrel_divd, Wrel_divs, out=np.zeros_like(Wrel_divs), where=Wrel_divs!=0)
    return Wrel, Wlower, Wupper    
    
def soil_evaporation_stage_two(th, th_s, th_fc, th_wp, th_dry, EsPot, EsAct, dz, dzsum, EvapZmin, EvapZmax, REW, EvapZ, fWrelExp, fevap, Wstage2, ToExtract, EvapTimeSteps):

    dims = th.shape
    nc, nr, nlat, nlon = dims[0], dims[1], dims[2], dims[3]
    thnew = np.copy(th)
    dz = dz[:,None,None,None] * np.ones((nr, nlat, nlon))[None,:,:,:]
    dzsum = dzsum[:,None,None,None] * np.ones((nr, nlat, nlon))[None,:,:,:]
    
    # Extract water
    cond11 = (ToExtract > 0)
    if np.any(cond11):
        Edt = ToExtract / EvapTimeSteps
        
        # Loop sub-daily time steps
        for jj in range(EvapTimeSteps):

            # Get water storage (mm) at start of stage 2 evaporation
            Wrel, Wlower, Wupper = relative_depletion(thnew, th_s, th_fc, th_wp, th_dry, dz[:,0,0,0], EvapZ, Wstage2, REW)

            # Check if need to expand evaporative layer
            cond111 = (cond11 & (EvapZmax > EvapZmin))
            Wcheck = (fWrelExp * ((EvapZmax - EvapZ) / (EvapZmax - EvapZmin)))                

            while np.any(cond111 & (Wrel < Wcheck) & (EvapZ < EvapZmax)):
                cond1111 = (cond111 & (Wrel < Wcheck) & (EvapZ < EvapZmax))

                # Expand evaporation layer by 1mm
                EvapZ[cond1111] += 0.001

                # Recalculate current water storage for new EvapZ
                Wrel, Wlower, Wupper = relative_depletion(thnew, th_s, th_fc, th_wp, th_dry, dz[:,0,0,0], EvapZ, Wstage2, REW)
                Wcheck = (fWrelExp * ((EvapZmax - EvapZ) / (EvapZmax - EvapZmin)))

            # Get stage 2 evaporation reduction coefficient
            Kr = ((np.exp(fevap * Wrel) - 1) / (np.exp(fevap) - 1))
            Kr = np.clip(Kr, None, 1)

            # Get water to extract (NB Edt is zero in cells which do not
            # need stage 2, so no need for index)
            ToExtractStg2 = (Kr * Edt)
            thnew, EsAct, ToExtract, ToExtractStg2 = extract_water(
                thnew, th_dry, dz, dzsum, EvapZmin, EsAct, ToExtract, ToExtractStg2)
            
    return thnew, EsAct
    
    
def adjust_potential_soil_evaporation_for_irrigation(surface_storage, EsPot, irr_method, irr_depth, prec, wet_surf):
    EsPotIrr = np.copy(EsPot)
    cond1 = ((irr_depth > 0) & (irr_method != 4))
    cond11 = (cond1 & ((prec > 0) | (surface_storage > 0)))
    EsPotIrr[cond11] = EsPot[cond11]
    cond12 = (cond1 & np.logical_not(cond11))
    EsPotIrr[cond12] = (EsPot * (wet_surf / 100))[cond12]  # TODO: more informative name for wet_surf
    return EsPotIrr
    
def adjust_potential_soil_evaporation_for_mulches(growing_season, surface_storage, EsPot, Mulches, fMulch, MulchPctGS, MulchPctOS):

    # NB if surface is flooded then there is no adjustment of potential soil
    # evaporation, regardless of mulches    
    EsPotMul = np.copy(EsPot)
    cond1 = (surface_storage < 0.000001)
    cond11 = (cond1 & (Mulches == 1))
    cond111 = (cond11 & growing_season)
    EsPotMul[cond111] = (EsPot * (1 - fMulch * (MulchPctGS / 100)))[cond111]
    cond112 = (cond11 & np.logical_not(growing_season))
    EsPotMul[cond112] = (EsPot * (1 - fMulch * (MulchPctOS / 100)))[cond112]
    return EsPotMul

def adjust_time_for_delayed_development(calendar_type, DAP, DelayedCDs, DelayedGDDs, GDDcum):
    # Adjust time for any delayed development
    if calendar_type == 1:
        tAdj = (DAP - DelayedCDs) # * growing_season_index
    elif calendar_type == 2:
        tAdj = (GDDcum - DelayedGDDs) # * growing_season_index
    return tAdj

def potential_soil_evaporation_rate(growing_season, tAdj, CC, CCadj, CCxAct, CCxW, Senescence, PrematSenes, et0, Kex, fwcc):

    dims = growing_season.shape
    nr, nlat, nlon = dims[0], dims[1], dims[2]
    # No canopy cover outside of growing season so potential soil
    # evaporation only depends on reference evapotranspiration
    EsPot = (Kex * et0)

    # Calculate maximum potential soil evaporation and potential soil
    # evaporation given current canopy size
    EsPotMax = (Kex * et0 * (1 - CCxW * (fwcc / 100)))
    EsPot[growing_season] = (Kex * (1 - CCadj) * et0)[growing_season]

    # Adjust potential soil evaporation for effects of withered canopy
    cond3 = (growing_season & (tAdj > Senescence) & (CCxAct > 0))
    mult = np.ones((nr, nlat, nlon))
    cond31 = (cond3 & (CC > (CCxAct / 2)))
    cond311 = (cond31 & (CC > CCxAct))
    mult[cond311] = 0
    cond312 = (cond31 & np.logical_not(cond311))
    mult_divd = (CCxAct - CC)
    mult_divs = (CCxAct / 2)
    mult[cond312] = np.divide(mult_divd, mult_divs, out=np.zeros_like(mult_divs), where=mult_divs!=0)[cond312]

    EsPot[cond3] = (EsPot * (1 - CCxAct * (fwcc / 100) * mult))[cond3]
    CCxActAdj = ((1.72 * CCxAct) + (CCxAct ** 2) - 0.3 * (CCxAct ** 3))
    EsPotMin = np.zeros((nr, nlat, nlon))
    EsPotMin[cond3] = (Kex * (1 - CCxActAdj) * et0)[cond3]
    EsPotMin = np.clip(EsPotMin, 0, None)

    # Line 85-89 of AOS_SoilEvaporation.m
    EsPot[cond3] = np.clip(EsPot, EsPotMin, EsPotMax)[cond3]

    cond4 = (growing_season & PrematSenes)
    EsPot[cond4] = np.clip(EsPot, None, EsPotMax)[cond4]
    return EsPot
    
def maximum_capillary_rise(ksat, aCR, bCR, zGW, zMid):
    """Function to compute maximum capillary rise for a given soil 
    profile
    
    Args:
      ksat : saturated hydraulic conductivity of the soil layer
      aCR  : ... of the soil layer
      bCR  : ... of the soil layer
      zGW  : depth of groundwater table
      z    : depth to midpoint of the soil layer

    """
    dims = ksat.shape
    nr, nlat, nlon = dims[0], dims[1], dims[2]
    MaxCR = np.zeros((nr, nlat, nlon))
    cond1 = ((ksat > 0) & (zGW > 0) & ((zGW - z) < 4))
    cond11 = (cond1 & (z >= zGW))
    MaxCR[cond11] = 99            
    cond12 = (cond1 & np.logical_not(cond11))
    MaxCR_log = np.log(zGW - z, out=np.zeros((MaxCR.shape)), where=cond12)
    MaxCR[cond12] = (np.exp((MaxCR_log - bCR) / aCR))[cond12]
    MaxCR = np.clip(MaxCR, None, 99)
    return MaxCR

def driving_force(th, th_fc_adj, th_wp, fshape_cr):
    """Function to calculate driving force of capillary rise
    
    Args:
      th        : water content of soil layer
      th_fc_adj : adjusted field capacity of soil layer
      th_wp     : wilting point of soil layer
      fshape_cr : ... of soil layer
    """
    dims = th.shape
    nr, nlat, nlon = dims[0], dims[1], dims[2]
    Df = np.ones((nr, nlat, nlon))
    cond = ((th >= th_wp) & (fshape_cr > 0))
    Df[cond] = (1 - (((th - th_wp) / (th_fc_adj - th_wp)) ** fshape_cr))[cond]
    Df = np.clip(Df, 0, 1)
    return Df

def relative_hydraulic_conductivity(th, th_fc, th_wp):
    dims = th.shape
    nr, nlat, nlon = dims[0], dims[1], dims[2]
    thThr = (th_wp + th_fc) / 2
    Krel = np.ones((nr, nlat, nlon))
    cond1 = (th < thThr)
    cond11 = (cond1 & ((th <= th_wp) | (thThr <= th_wp)))
    Krel[cond11] = 0
    cond12 = (cond1 & np.logical_not(cond11))
    Krel[cond12] = ((th - th_wp) / (thThr - th_wp))[cond12]
    return Krel

def store_water_from_capillary_rise(th, th_fc, th_fc_adj, th_wp, fshape_cr, MaxCR, flux_out, zGW, zBot, dz):
    dims = th.shape
    nc, nlat, nlon = dims[0], dims[1], dims[2]
    thnew = np.copy(th)

    cond1 = ((np.round(MaxCR * 1000) > 0) & (np.round(flux_out * 1000) == 0))

    # calculate driving force
    Df = driving_force(th, th_fc_adj, th_wp, fshape_cr)

    # Calculate relative hydraulic conductivity
    Krel = relative_hydraulic_conductivity(th, th_fc, th_wp)

    # Check if room is available to store water from capillary rise
    dth = np.zeros((nr, nlat, nlon))
    dth[cond1] = (th_fc_adj - th)[cond1]                
    dth = np.round((dth * 1000) / 1000)

    # Store water if room is available
    dthMax = np.zeros((nr, nlat, nlon))
    CRcomp = np.zeros((nr, nlat, nlon))
    cond15 = (cond1 & (dth > 0) & ((zBot - dz / 2) < zGW))

    dthMax[cond35] = (Krel * Df * MaxCR / (1000 * dz))[cond35]
    cond151 = (cond15 & (dth >= dthMax))
    thnew[cond151] += dthMax[cond151]
    CRcomp[cond151] = (dthMax * 1000 * dz)[cond151]
    MaxCR[cond151] = 0

    cond152 = (cond35 & np.logical_not(cond151))
    thnew[cond152] = th_fc_adj[cond152]
    CRcomp[cond152] = (dth * 1000 * dz)[cond152]
    MaxCR[cond152] = ((Krel * MaxCR) - CRcomp)[cond152]    
    return thnew, CRcomp, MaxCR

def capillary_rise(th, th_fc, th_fc_adj, th_wp, fshape_cr, aCR, bCR, MaxCR, flux_out, zGW, dz):
    
    dims = th.shape
    nc, nr, nlat, nlon = dims[0], dims[1], dims[2], dims[3]
    zBot = np.sum(dz)
    thnew = np.copy(th)
    compi = nc - 1
    WCr = np.zeros((nr, nlat, nlon))
    while ((np.any(np.round(MaxCR * 1000) > 0)) & (np.any(np.round(flux_out[compi,:] * 1000) == 0)) & (compi >= 0)):

        cond1 = ((np.round(MaxCR * 1000) > 0) & (np.round(flux_out[compi,:] * 1000) == 0))
        thnew_comp, CRcomp, MaxCR = store_water_from_capillary_rise(
            th[compi,:], th_fc[compi,:], th_fc_adj[compi,:], th_wp[compi,:],
            fshape_cr, MaxCR, flux_out[compi,:], zGW, zBot, dz[compi])

        thnew[compi,:][cond1] = thnew_comp[cond1]
        WCr[cond1] += CRcomp[cond1]

        # Update bottom elevation of compartment and compartment counter
        zBot -= dz[compi]
        compi -= 1

        # Update restriction on maximum capillary rise
        if compi >= 0:
            zBotMid = zBot - (dz[compi] / 2)
            LimCR = maximum_capillary_rise(ksat[compi,:], aCR[compi,:], bCR[compi,:], zGW, zBotMid)
            cond11 = (cond1 & (MaxCR > LimCR))
            MaxCR[cond11] = LimCR[cond11]

    return thnew, WCr
    
def infiltration(to_store, flux_out, th, th_s, th_fc, th_fc_adj, tau, ksat, dz):
    """Function to infiltrate incoming water into the soil column"""
    
    thnew = np.copy(th)
    dims = th.shape
    nc, nr, nlat, nlon = dims[0], dims[1], dims[2], dims[3]

    # Initialize counters
    comp = 0
    runoff = np.zeros((nr, nlat, nlon))
    cond3_ini = (to_store > 0)

    while (np.any(to_store > 0) & (comp < nc)):

        # Update condition
        cond3 = (to_store > 0)

        # Calculate saturated drainage ability, drainage factor and
        # required drainage ability
        dthdtS = tau[comp,:] * (th_s[comp,:] - th_fc[comp,:])
        factor = ksat[comp,:] / (dthdtS * 1000 * dz[comp])
        dthdt0 = to_store / (1000 * dz[comp])

        # Initialize water content for current layer
        theta0 = np.zeros((nr, nlat, nlon))
        
        # Check drainage ability
        cond31 = (cond3 & (dthdt0 < dthdtS))

        # Calculate water content needed to meet drainage dthdt0
        cond311 = (cond31 & (dthdt0 <= 0))
        theta0[cond311] = th_fc_adj[comp,:][cond311]
        cond312 = (cond31 & np.logical_not(cond311))
        A = (1 + ((dthdt0 * (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)) / (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]))))
        theta0[cond312] = (th_fc[comp,:] + np.log(A))[cond312]

        # Limit thX to between saturation and field capacity
        cond313 = (cond31 & (theta0 > th_s[comp,:]))
        theta0[cond313] = th_s[comp,:][cond313]
        cond314 = (cond31 & np.logical_not(cond313) & (theta0 < th_fc_adj[comp,:]))
        theta0[cond314] = th_fc_adj[comp,:][cond314]
        dthdt0[cond314] = 0

        # Limit water content and drainage to saturation
        cond32 = (cond3 & np.logical_not(cond31))
        theta0[cond32] = th_s[comp,:][cond32]
        dthdt0[cond32] = dthdtS[cond32]

        # Calculate maximum water flow through compartment and total
        # drainage
        drainmax = factor * dthdt0 * 1000 * dz[comp]
        drainage = drainmax + flux_out[comp,:]

        # Limit drainage to saturated hydraulic conductivity
        cond33 = (cond3 & (drainage > ksat[comp,:]))
        drainmax[cond33] = (ksat[comp,:] - flux_out[comp,:])[cond33]

        # Line 117 of AOS_Infiltration.m
        # Calculate difference between threshold and current water contents
        diff = theta0 - th[comp,:]

        cond34 = (cond3 & (diff > 0))
        thnew[comp,:][cond34] += (to_store / (1000 * dz[comp]))[cond34]

        # Remaining water that can infiltrate to compartments below
        cond341 = (cond34 & (thnew[comp,:] > theta0))
        to_store[cond341] = ((thnew[comp,:] - theta0) * 1000 * dz[comp])[cond341]
        thnew[comp,:][cond341] = theta0[cond341]

        # Otherwise all infiltrating water has been stored
        cond342 = (cond34 & np.logical_not(cond341))
        to_store[cond342] = 0

        # Update outflow from current compartment (drainage + infiltration
        # flows)
        flux_out[comp,:][cond3] += to_store[cond3]

        # Calculate backup of water into compartments above, and update
        # water to store
        excess = np.clip((to_store - drainmax), 0, None)
        to_store -= excess

        # Redistribute excess to compartments above
        precomp = comp + 1
        while (np.any(cond3 & (excess > 0)) & (precomp != 0)):

            # Update condition
            cond35 = (cond3 & (excess > 0))

            # Keep storing in compartments above until soil surface is
            # reached
            precomp -= 1

            # Update outflow from compartment
            flux_out[precomp,:][cond35] = (flux_out[precomp,:] - excess)[cond35]

            # Update water content and limit to saturation
            thnew[precomp,:][cond35] += (excess / (1000 * dz[precomp]))[cond35]
            cond351 = (cond35 & (thnew[precomp,:] > th_s[precomp,:]))
            excess[cond351] = ((thnew[precomp,:] - th_s[precomp,:]) * 1000 * dz[precomp])[cond351]
            thnew[precomp,:][cond351] = th_s[precomp,:][cond351]
            cond352 = (cond35 & np.logical_not(cond351))
            excess[cond352] = 0

        # Any leftover water not stored becomes runoff
        cond36 = (cond3 & (excess > 0))
        runoff[cond36] += excess[cond36]

        # update comp
        comp += 1

    # Infiltration left to store after bottom compartment becomes deep
    # percolation (mm)
    deep_perc = np.copy(to_store)

    # Otherwise if to_store equals zero there is no infiltration
    cond4 = np.logical_not(cond3_ini)
    deep_perc[cond4] = 0
    runoff[cond4] = 0
    return runoff, flux_out, deep_perc, thnew
        
def update_surface_storage(bunds, zbund, ksat, infl, surface_storage0):
    """Function to update surface storage, infiltration storage and 
    runoff
    """
    dims = ksat.shape
    nc, nr, nlon, nlat = dims[0], dims[1], dims[2], dims[3]
    to_store = np.zeros((nr, nlat, nlon))
    surface_storage = np.copy(surface_storage0)
    runoff = np.zeros((nr, nlat, nlon))
    infl_tot = infl + surface_storage0
    
    # Infiltration limited by saturated hydraulic conductivity of surface
    # soil layer; additional water ponds on surface
    cond1 = (bunds == 1) & (zbund > 0.001)
    cond11 = (cond1 & (infl_tot > 0))
    cond111 = (cond11 & (infl_tot > ksat[0,:]))
    to_store[cond111] = ksat[0,:][cond111]
    surface_storage[cond111] = (infl_tot - ksat[0,:])[cond111]

    # Otherwise all water infiltrates and surface storage becomes zero
    cond112 = (cond11 & np.logical_not(cond111))
    to_store[cond112] = infl_tot[cond112]
    surface_storage[cond112] = 0

    # Calculate additional runoff if water overtops bunds
    cond113 = (cond11 & (surface_storage > (zbund * 1000)))
    runoff[cond113] = (surface_storage - (zbund * 1000))[cond113]
    surface_storage[cond113] = (zbund * 1000)[cond113]

    # Otherwise excess water does not overtop bunds and there is no runoff
    cond114 = (cond11 & np.logical_not(cond113))
    runoff[cond114] = 0

    # If total infiltration is zero then there is no storage or runoff
    cond12 = (cond1 & np.logical_not(cond11))
    to_store[cond12] = 0
    runoff[cond12] = 0

    # If there are no bunds then infiltration is divided between runoff
    # and infiltration according to saturated conductivity of surface
    # layer
    cond2 = (bunds == 0)
    cond21 = (cond2 & (infl > ksat[0,:]))
    to_store[cond21] = ksat[0,:][cond21]
    runoff[cond21] = (infl - ksat[0,:])[cond21]
    cond22 = (cond2 & np.logical_not(cond21))
    to_store[cond22] = infl[cond22]
    runoff[cond22] = 0
    return surface_storage, runoff, to_store
    
def rainfall_partition(P, th, th_fc, th_wp, zcn, adjcn, cn, cnbot, cntop, dz, dzsum, bunds, zbund):

    dims = th_fc.shape
    nc, nr, nlon, nlat = dims[0], dims[1], dims[2], dims[3]
    
    # Add rotation dimension to precipitation
    P = P[None,:,:] * np.ones((nr))[:,None,None]

    arr_ones = np.ones((nr, nlon, nlat))[None,:,:,:]
    dz = dz[:,None,None,None] * arr_ones
    dzsum = dzsum[:,None,None,None] * arr_ones
    zcn = zcn[None,:,:,:] * np.ones((nc))[:,None,None,None]

    Runoff = np.zeros((nr, nlon, nlat))
    Infl = np.zeros((nr, nlon, nlat))
    
    cond1 = ((bunds == 0) | (zbund < 0.001))

    # Adjust curve number
    cond11 = (cond1 & (adjcn == 1))

    # Check which compartment cover depth of top soil used to adjust
    # curve number
    comp_sto = (np.round(dzsum * 1000) <= np.round(zcn * 1000))
    cond111 = np.all((comp_sto == False), axis=0)
    cond111 = np.broadcast_to(cond111, comp_sto.shape)
    comp_sto[cond111] = True

    # Calulcate weighting factors by compartment
    dzsum[dzsum > zcn] = zcn[dzsum > zcn]
    wx = (1.016 * (1 - np.exp(-4.16 * (dzsum / zcn))))

    # xx is wx for the overlying layer, with the top layer equal to zero
    xx = np.concatenate((np.zeros((1, nr, nlat, nlon)), wx[:-1,:]), axis=0)
    wrel = np.clip((wx - xx), 0, 1)

    # Calculate relative wetness of topsoil
    thnew = np.maximum(th_wp, th)

    # Multiply by comp_sto to ensure that compartments not used for
    # curve number adjustment are set to zero
    wet_top_comp = wrel * ((th - th_wp) / (th_fc - th_wp)) * comp_sto
    wet_top = np.sum(wet_top_comp, axis=0)
    wet_top = np.clip(wet_top, 0, 1)

    # Make adjustment to curve number
    # cn = np.copy(cn)
    cn[cond11] = np.round(cnbot + (cntop - cnbot) * wet_top)[cond11]

    # Partition rainfall into runoff and infiltration
    S = (25400. / cn) - 254
    term = (P - ((5. / 100) * S))

    cond12 = (cond1 & (term <= 0))
    Runoff[cond12] = 0
    Infl[cond12] = P[cond12]

    cond13 = (cond1 & np.logical_not(cond12))
    Runoff[cond13] = ((term ** 2) / (P + (1 - (5. / 100)) * S))[cond13]
    Infl[cond13] = (P - Runoff)[cond13]

    # If bunds are present then there is no runoff
    cond2 = np.logical_not(cond1)
    Runoff[cond2] = 0
    Infl[cond2] = P[cond2]
    return Runoff,Infl

def compute_dthdt(th, th_s, th_fc, th_fc_adj, tau):
    dims = th.shape
    dthdt = np.zeros((dims[0], dims[1], dims[2]))
    cond1 = th <= th_fc_adj
    dthdt[cond1] = 0
    cond2 = np.logical_not(cond1) & (th >= th_s)
    dthdt[cond2] = (tau * (th_s - th_fc))[cond2]
    cond3 = np.logical_not(cond1 | cond2)
    dthdt[cond3] = (tau * (th_s - th_fc) * ((np.exp(th - th_fc) - 1) / (np.exp(th_s - th_fc) - 1)))[cond3]
    cond4 = ((cond2 | cond3) & ((th - dthdt) < th_fc_adj))
    dthdt[cond4] = (th - th_fc_adj)[cond4]
    return dthdt
        
def drainage(th, th_s, th_fc, th_fc_adj, ksat, tau, dz, dzsum):

    dims = th.shape    
    thnew = np.copy(th)
    drainsum = np.zeros((dims[1], dims[2], dims[3]))
    FluxOut = np.zeros((dims[0], dims[1], dims[2], dims[3]))
    for comp in range(dims[0]):

        # Calculate drainage ability of compartment ii
        dthdt = compute_dthdt(th[comp,:], th_s[comp,:], th_fc[comp,:], th_fc_adj[comp,:], tau[comp,:])

        # Drainage from compartment ii (mm) (Line 41 in AOS_Drainage.m)
        draincomp = dthdt * dz[comp] * 1000

        # Check drainage ability of compartment ii against cumulative
        # drainage from compartments above (Lines 45-52 in AOS_Drainage.m)
        excess = np.zeros((dims[1], dims[2], dims[3]))
        prethick = dzsum[comp] - dz[comp]
        drainmax = dthdt * 1000 * prethick
        drainability = (drainsum <= drainmax)

        # Drain compartment
        cond5 = drainability
        thnew[comp,:][cond5] = (th[comp,:] - dthdt)[cond5]

        # Update cumulative drainage (mm), restrict to saturated hydraulic
        # conductivity and adjust excess drainage flow
        drainsum[cond5] += draincomp[cond5]
        cond51 = (cond5 & (drainsum > ksat[comp,:]))
        excess[cond51] += (drainsum - ksat[comp,:])[cond51]
        drainsum[cond51] = ksat[comp,:][cond51]

        # Calculate value of theta (thX) needed to provide a drainage
        # ability equal to cumulative drainage (Lines 70-85 in AOS_Drainage.m)
        cond6 = np.logical_not(drainability)
        dthdt[cond6] = np.divide(drainsum, 1000 * prethick, out=np.zeros_like(drainsum), where=prethick!=0)[cond6]
        thX = np.zeros((dims[1], dims[2], dims[3]))
        cond61 = (cond6 & (dthdt <= 0))
        thX[cond61] = th_fc_adj[comp,:][cond61]
        cond62 = (cond6 & np.logical_not(cond61) & (tau[comp,:] > 0))
        A = (1 + ((dthdt * (np.exp(th_s[comp,:] - th_fc[comp,:]) - 1)) / (tau[comp,:] * (th_s[comp,:] - th_fc[comp,:]))))
        thX[cond62] = (th_fc_adj[comp,:] + np.log(A))[cond62]
        thX[cond62] = np.clip(thX, th_fc_adj[comp,:], None)[cond62]
        cond63 = (cond6 & np.logical_not(cond61 | cond62))
        thX[cond63] = (th_s[comp,:] + 0.01)[cond63]

        # Check thX against hydraulic properties of current soil layer

        # Increase compartment ii water content with cumulative drainage
        cond64 = (cond6 & (thX <= th_s[comp,:]))
        thnew[comp,:][cond64] = (th[comp,:] + (drainsum / (1000 * dz[comp])))[cond64]

        # Cumulative drainage is the drainage difference between theta_x and
        # new theta plus drainage ability at theta_x
        cond641 = (cond64 & (thnew[comp,:] > thX))
        drainsum[cond641] = ((thnew[comp,:] - thX) * 1000 * dz[comp])[cond641]

        # Calculate drainage ability for thX
        dthdt = compute_dthdt(thX, th_s[comp,:], th_fc[comp,:], th_fc_adj[comp,:], tau[comp,:])

        # Update cumulative drainage (mm), restrict to saturated hydraulic
        # conductivity and adjust excess drainage flow
        drainsum[cond641] += (dthdt * 1000 * dz[comp])[cond641]
        cond6415 = (cond641 & (drainsum > ksat[comp,:]))
        excess[cond6415] += (drainsum - ksat[comp,:])[cond6415]
        drainsum[cond6415] = ksat[comp,:][cond6415]

        # Update water content
        thnew[comp,:][cond641] = (thX - dthdt)[cond641]

        # Calculate drainage ability for updated water content
        cond642 = (cond64 & np.logical_not(cond641) & (thnew[comp,:] > th_fc_adj[comp,:]))
        dthdt = compute_dthdt(thnew[comp,:], th_s[comp,:], th_fc[comp,:], th_fc_adj[comp,:], tau[comp,:])

        # Update water content
        thnew[comp,:][cond642] = (thnew[comp,:] - dthdt)[cond642]

        # Update cumulative drainage (mm), restrict to saturated hydraulic
        # conductivity and adjust excess drainage flow
        drainsum[cond642] = (dthdt * 1000 * dz[comp])[cond642]
        cond6425 = (cond642 & (drainsum > ksat[comp,:]))
        excess[cond6425] += (drainsum - ksat[comp,:])[cond6425]
        drainsum[cond6425] = ksat[comp,:][cond6425]

        # Otherwise, drainage is zero
        cond643 = (cond64 & np.logical_not(cond641 | cond642))
        drainsum[cond643] = 0

        # Increase water content in compartment ii with cumulative
        # drainage from above
        cond65 = (cond6 & np.logical_not(cond64) & (thX > th_s[comp,:]))
        thnew[comp,:][cond65] = (th[comp,:] + (drainsum / (1000 * dz[comp])))[cond65]

        # Check new water content against hydraulic properties of soil layer
        # Lines 166-198
        cond651 = (cond65 & (thnew[comp,:] <= th_s[comp,:]))

        # Calculate new drainage ability
        cond6511 = (cond651 & (thnew[comp,:] > th_fc_adj[comp,:]))
        dthdt = compute_dthdt(thnew[comp,:], th_s[comp,:], th_fc[comp,:], th_fc_adj[comp,:], tau[comp,:])

        # Update water content
        thnew[comp,:][cond6511] -= (dthdt)[cond6511]

        # Update cumulative drainage (mm), restrict to saturated hydraulic
        # conductivity and adjust excess drainage flow
        drainsum[cond6511] = (dthdt * 1000 * dz[comp])[cond6511]
        cond65115 = (cond6511 & (drainsum > ksat[comp,:]))
        excess[cond65115] += (drainsum - ksat[comp,:])[cond65115]
        drainsum[cond65115] = ksat[comp,:][cond65115]

        cond6512 = (cond651 & (np.logical_not(cond6511)))
        drainsum[cond6512] = 0

        # Calculate excess drainage above saturation
        cond652 = (cond65 & np.logical_not(cond651) & (thnew[comp,:] > th_s[comp,:]))
        excess[cond652] = ((thnew[comp,:] - th_s[comp,:]) * 1000 * dz[comp])[cond652]

        # Calculate drainage ability for updated water content
        dthdt = compute_dthdt(thnew[comp,:], th_s[comp,:], th_fc[comp,:], th_fc_adj[comp,:], tau[comp,:])

        # Update water content
        thnew[comp,:][cond652] = (th_s[comp,:] - dthdt)[cond652]

        # Update drainage, maximum drainage, excess drainage
        draincomp[cond652] = (dthdt * 1000 * dz[comp])[cond652]
        drainmax[cond652] = (dthdt * 1000 * prethick)[cond652]
        drainmax[cond652] = np.clip(drainmax, None, excess)[cond652]
        excess[cond652] -= drainmax[cond652]

        # Update cumulative drainage (mm), restrict to saturated hydraulic
        # conductivity and adjust excess drainage flow
        drainsum[cond652] = (draincomp + drainmax)[cond652]
        cond6525 = (cond652 & (drainsum > ksat[comp,:]))
        excess[cond6525] += (drainsum - ksat[comp,:])[cond6525]
        drainsum[cond6525] = ksat[comp,:][cond6525]

        # Store output flux from compartment ii
        FluxOut[comp,:] = np.copy(drainsum)

        # Redistribute excess in compartment above
        precomp = comp + 1
        while (np.any(excess > 0)) & (precomp != 0):

            # Include condition here so that it is updated
            cond7 = (excess > 0)

            # Update compartment counter
            precomp -= 1

            # Update flux from compartment
            if (precomp < comp):
                FluxOut[precomp,:][cond7] -= excess[cond7]

            # Increase water content to store excess
            thnew[precomp,:][cond7] += (excess / (1000 * dz[precomp]))[cond7]

            # Limit water content to saturation and adjust excess counter
            cond71 = (cond7 & (thnew[precomp,:] > th_s[precomp,:]))
            excess[cond71] = ((thnew[precomp,:] - th_s[precomp,:]) * 1000 * dz[precomp])[cond71]
            thnew[precomp,:][cond71] = th_s[precomp,:][cond71]

            cond72 = (cond7 & np.logical_not(cond71))
            excess[cond72] = 0

    # Update state variables        
    DeepPerc = np.copy(drainsum)
    th = np.copy(thnew)
    return th, DeepPerc, FluxOut

def root_zone_water(th, th_s, th_fc, th_wp, th_dry, dz, zmin, zroot, aer):

    dims = th.shape
    arr_ones = np.ones((dims[1], dims[2], dims[3]))
    dz = dz[:,None,None,None] * arr_ones
    dzsum = np.cumsum(dz, axis=0)

    # Calculate root zone water content and available water
    rootdepth = np.maximum(zmin, zroot)
    rootdepth = np.round(rootdepth * 100) / 100
    comp_sto = (np.round((dzsum - dz) * 1000) < np.round(rootdepth * 1000))

    # Fraction of compartment covered by root zone (zero in compartments
    # NOT covered by the root zone)
    factor = 1 - ((dzsum - rootdepth) / dz)
    factor = np.clip(factor, 0, 1)
    factor[np.logical_not(comp_sto)] = 0

    # Water storages in root zone (mm) - initially compute value in each
    # compartment, then sum to get overall root zone storages
    Wr_comp = factor * 1000 * th * dz
    WrS_comp = factor * 1000 * th_s * dz
    WrFC_comp = factor * 1000 * th_fc * dz
    WrWP_comp = factor * 1000 * th_wp * dz
    WrDry_comp = factor * 1000 * th_dry * dz

    # Water storage in root zone at aeration stress threshold (mm)
    WrAer_comp = factor * 1000 * (th_s - (aer / 100)) * dz

    Wr = np.sum(Wr_comp, axis=0)
    Wr[Wr < 0] = 0
    WrS = np.sum(WrS_comp, axis=0)
    WrFC = np.sum(WrFC_comp, axis=0)
    WrWP = np.sum(WrWP_comp, axis=0)
    WrDry = np.sum(WrDry_comp, axis=0)
    WrAer = np.sum(WrAer_comp, axis=0)

    # Convert depths to m3/m3
    thRZ = {}
    thRZ['Act'] = np.divide(Wr, rootdepth * 1000, out=np.zeros_like(Wr), where=rootdepth!=0)
    thRZ['Sat'] = np.divide(WrS, rootdepth * 1000, out=np.zeros_like(WrS), where=rootdepth!=0)
    thRZ['Fc']  = np.divide(WrFC, rootdepth * 1000, out=np.zeros_like(WrFC), where=rootdepth!=0)
    thRZ['Wp']  = np.divide(WrWP, rootdepth * 1000, out=np.zeros_like(WrWP), where=rootdepth!=0)
    thRZ['Dry'] = np.divide(WrDry, rootdepth * 1000, out=np.zeros_like(WrDry), where=rootdepth!=0)
    thRZ['Aer'] = np.divide(WrAer, rootdepth * 1000, out=np.zeros_like(WrAer), where=rootdepth!=0)

    # Calculate total available water and root zone depletion
    TAW = np.clip((WrFC - WrWP), 0, None)
    Dr = np.clip((WrFC - Wr), 0, None)
    Wr = np.copy(Wr)
    return thRZ, TAW, Dr, Wr
        
