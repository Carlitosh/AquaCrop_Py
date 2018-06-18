#!/usr/bin/env python
# -*- coding: utf-8 -*-

# AquaCrop crop growth model

import os
import shutil
import sys
import math
import gc
import numpy as np

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
    
