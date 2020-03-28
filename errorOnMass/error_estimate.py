#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from astropy.coordinates import SkyCoord
import vector_astrometry as va
import astrometry as a
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.io.fits import getdata
import leigh as l


# separation at a given julian year
def sepatt(jyear):
    time = Time(jyear, format='jyear')
    l = va.T_pos((time, t0, lcoord), lpmra, lpmdec, lplx, lvrad)
    s = va.T_pos((time, t0, scoord), spmra, spmdec, splx, svrad)
    return l.separation(s).arcsec * 1000.


# Returns mass given distance to lens, true separation and deflection in pc, mas, mas respectively
def M(D, s, d):
    return D* np.power((np.sqrt(d) * np.sqrt(d+s) / 90.25), 2)




# deflection calculation
def defl(sep):
    M = 0.61
    Dl = 4.63
    Oe = 90.25 * np.sqrt(M/Dl)
    u = sep/Oe
    return 0.5 * (np.sqrt(np.power(u,2) + 4) - u ) * Oe



# tangent plane coords of lens at given times
def tpcoordatT__l(time, raz, decz):
    l = va.T_pos((time, t0, lcoord), lpmra, lpmdec, lplx, lvrad)
    lxi, lxn, adz = a.s2tp(l.ra.deg, l.dec.deg, raz=raz, decz=decz, unit="deg")
    return lxi, lxn, adz
    
    
# tangent plane coords of source at given times
def tpcoordatT__s(time, raz, decz):
    s = va.T_pos((time, t0, scoord), spmra, spmdec, splx, svrad)
    sxi, sxn, adz = a.s2tp(s.ra.deg, s.dec.deg, raz=raz, decz=decz, unit="deg")
    return sxi, sxn, adz
    
    

if __name__=="__main__":
    ax = plt.axes()

    # get gost data
    gost = getdata('gost_LAWD37.fits', -1, view=np.recarray)

    #########################
    ### For 'Candidate 4' ###
    #########################
    t0 = Time(2015.0, format='jyear')
    
    lcoord = SkyCoord(ra=176.4549073, dec=-64.84295714, unit='deg', frame='icrs')
    lpmra = 2.66203572627 # arcsec/yr
    lpmdec = -0.34518255501 # arcsec/yr
    lplx = 0.2157823335 # arcsec
    lvrad = 0.0 # km/s
    
    scoord = SkyCoord(ra=176.463605, dec=-64.8432977, unit='deg', frame='icrs')
    spmra = -14.3 * 1E-3  # arcsec/yr
    spmdec = -2.0 * 1E-3  # arcsec/yr
    splx = 2E-6 # arcsec, some tiny value
    svrad = 0.0 # km/s
    
    
    # times of Gaia transits
    t_gost = Time(gost['TIME_UTC_NAME'], format='isot')
    

    # tangent plane coords at Gaia transit times
    glxi, glxn, adz = tpcoordatT__l(t_gost, None, None)
    gsxi, gsxn, _ = tpcoordatT__s(t_gost, adz[0], adz[1])
    
    # x and y separations
    x_sep = (glxi - gsxi)*1000.
    y_sep = (glxn - gsxn)*1000.
    
    # total true separation
    s = np.hypot(x_sep, y_sep)
    # deflection given true separation
    d = defl(s)
    
    # Gaia estimated error on single epoch distance
    e_d = 0.24
    e_s = 0.0
    e_D = 0.0
    
    # use only given sigma measurements
    these = d/e_d > 2.0
    s = s[these]
    d = d[these]
    
    
    runs = 1E6
    M_op, e_M_op = [], []
    for _s, _d in zip(s, d):
        _s = np.random.normal(loc=_s, scale=e_s, size=int(runs))
        _d = np.random.normal(loc=_d, scale=e_d, size=int(runs))
        _D = np.random.normal(loc=4.63, scale=e_D, size=int(runs))
        M_sample = M(_D, _s, _d)
        M_op += [np.nanmean(M_sample)]
        e_M_op += [np.nanstd(M_sample)]
    
    meas_M, e_meas_M = l.weightedMean(M_op, e_M_op)
    
    print("input mass: 0.61 Msun")
    print("recovered mass: %4.2f +- %4.2f Msun" % (meas_M, e_meas_M))
    print("%6.2f percent error" % (e_meas_M/meas_M * 100.))
    
    
    
    
