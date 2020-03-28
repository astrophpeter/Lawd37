#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy import sin, cos, radians
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel
from astropy import constants as const


def vector(x, y, z):
    return np.array([x, y, z]).T

def vecmag(vec):
    return np.sqrt(np.sum(np.power(vec, 2), axis=1))

def vecnorm(vec):
    _vecmag = vecmag(vec)
    _vecmag = np.repeat(_vecmag, vec.shape[1]).reshape((-1,3))
    return vec/_vecmag


def plx2au(plx):
    return 648000. / (np.pi * plx)


def R(time, vel=False):
    # R is the position vector of Earth
    # requires an astropy Time object
    # returns R, a three component array describing the position vector of Earth
    # in au
    _R, _Rdot = get_body_barycentric_posvel('earth', time, ephemeris='jpl')
    c = vector(_R.x.to(u.AU) / u.AU, _R.y.to(u.AU) / u.AU, _R.z.to(u.AU) / u.AU)
    if vel:
        d = vector(_Rdot.x.to(u.AU/u.d) / (u.AU/u.d),
                   _Rdot.y.to(u.AU/u.d) / (u.AU/u.d),
                   _Rdot.z.to(u.AU/u.d) / (u.AU/u.d))
        return c, d
    else:
        return c



def s(coord):
    # Return the unit vector of coordinates
    # Requires an astropy coordinate object(s)
    a, d = coord.ra.rad, coord.dec.rad
    return vector(
                   cos(a) * cos(d),
                   sin(a) * cos(d),
                   sin(d)
                 )


def m(coord, pmra, pmdec, plx, vr):
    # Return the space motion vector of a star
    # Requirements:
    #   - astropy coordinate object
    #   - proper motion in ra and dec (arcsec per year)
    #   - parallax (arcsec)
    #   - radial velocity (km/s)
    # Returns astropy cartesian representation of space motion vector
    # see eqns.
    #   - 11.2
    #   - 12.36
    # in 'Spherical Astronomy' by R.M. Green
    a, d = coord.ra.rad, coord.dec.rad

    # convert radial velocity and proper motion to AU/yr
    vr /= 4.74
    vt_ra = pmra / plx
    vt_dec = pmdec / plx

    # Tangential velocity vector
    # removed cos(d) terms on the vt_ra sections, I'm sure they are unnecessary
    # when pmra = pmra*cos(d)
    # I can't think why Green included the sin 1" terms, I've omitted them and
    # it works absolutely fine. sin 1" is tiny, it must have something to do
    # with the coordinates he expects or the operations he performs with this
    # vector and how I do them differently?
    Vtan = vector(
                  - vt_ra * sin(a) - vt_dec * cos(a) * sin(d),
                  vt_ra * cos(a) - vt_dec * sin(a) * sin(d),
                  vt_dec * cos(d)
                 )

    # Radial velocity vector
    Vrad = vr * s(coord)
    return Vtan + Vrad


def T_pos(t, t0, coords, pmra, pmdec, plx, vrad=0.0):
    # Return equatorial coordinates at times 't'
    # Requirements:
    #   - astropy time objects 't' and 't0', the output time(s) and start
    #     time respectively
    #   - astropy coordinate object describing the position of the target at
    #     time t0
    #   - proper motion of the target
    #   - parallax of the target
    #   - (optional) radial velocity of the target
    # The units of plx, pmra and pmdec should be arcsec(/year). The units
    # of vrad should be km/s

    # time difference (julian years)
    dt = t.jyear - t0.jyear
    if isinstance(dt, float):
        count = 1
    else:
        count = dt.shape[0]

    # start position vector of star
    _X0 = s(coords) * plx2au(plx)

    # space motion vector of star (in AU/yr)
    _m = m(coords, pmra, pmdec, plx, vrad)
    _m = np.tile(_m, count).reshape((-1,3)).T
    mdt = (_m*dt).T

    # position and velocity vectors of earth (in AU, AU/day)
    _R = R(t, vel=False)

    # position vector of star at time t
    _X1 = _X0 + mdt - _R

    # set astropy sky coordinates of position uncorrected for aberration, we
    # want the distance
    _s1 = SkyCoord(_X1[:,0], _X1[:,1], _X1[:,2], frame='icrs', unit='AU', representation='cartesian')
    _s1.representation = 'spherical'
    
    
    ## |||| not bothering with aberration ||||
    ## VVVV                               VVVV
    '''
    # c in AU per day
    c_aud = const.c.to(u.AU/u.d) / u.AU * u.d
    # position unit vector of star at time t
    s1 = vecnorm(_X1)
    # position unit vector of star corrected for aberration
    _s2 = vecnorm(s1 + _Rdot / c_aud)

    # set astropy sky coordinates (with false units for now since we dont care
    # about distance for the aberration correction)
    _s2 = SkyCoord(_s2[:,0], _s2[:,1], _s2[:,2], frame='icrs', representation='cartesian')
    _s2.representation = 'spherical'
    # give the uncorrected distance to the corrected position
    #_s2.distance = _s1.distance
    '''
    
    if count == 1:
        return _s1[0]
    else:
        return _s1



if __name__=="__main__":
    #########################
    ### For 'Candidate 4' ###
    #########################
    coord = SkyCoord(ra=176.45490, dec=-64.842955, unit='deg', frame='icrs')
    t0 = Time(2015.0, format='jyear')
    pmra = 2.66203572627
    pmdec = -0.34518255501
    plx = 0.2157823335
    vrad = 0.0

    tx = Time([2458465.5], format='jd')
    a,b = R(tx)
    #print a[0]
    #print b[0]

    ricky = np.genfromtxt('Ricky_lens.dat', skip_header=1)
    t = Time(ricky[:,0], format='jd')

    #t = Time(2019.9, format='jyear')

    c1 = T_pos((t, t0, coord), pmra, pmdec, plx, vrad)


    Ricky_coord = SkyCoord(ra=ricky[0,1], dec=ricky[0,4], unit='deg', frame='icrs')
    #print 'new - ricky:', c1.separation(Ricky_coord).arcsec[0] * 1000.

    #print 'new - ricky (acosd):', (c1.ra.deg[0] - Ricky_coord.ra.deg) * np.cos(coord.dec.rad) * 3E6
    #print 'new - ricky (d):', (c1.dec.deg[0] - Ricky_coord.dec.deg)*3E6


    from skyfield.api import load, Star, Angle
    lens = Star(ra=Angle(degrees=coord.ra.deg),
                dec=Angle(degrees=coord.dec.deg),
                ra_mas_per_year=pmra*1000.,
                dec_mas_per_year=pmdec*1000.,
                parallax_mas=plx*1000.)
    ts = load.timescale()
    ts1 = ts.utc(t.jyear[0]-15.)
    planets = load('de421.bsp')
    earth = planets['earth']
    astrometric = earth.at(ts1).observe(lens)
    rai, deci, distancesa = astrometric.radec()
    skyfield_coord = SkyCoord(ra=rai._degrees, dec=deci._degrees, unit='deg', frame='icrs')
    #print 'skyfield - ricky:', skyfield_coord.separation(Ricky_coord).arcsec * 1000.
    #print 'skyfield - new:', skyfield_coord.separation(c1[0]).arcsec * 1000.
