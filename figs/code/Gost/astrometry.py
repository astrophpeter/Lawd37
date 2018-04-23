#!/usr/bin/env python

import numpy as np
from numpy import sin, cos, radians
from astropy.time import Time
from astropy.coordinates import SkyCoord
from vector_astrometry import R


def cos2(x):
    return np.power(cos(x), 2)


def sin2(x):
    return np.power(sin(x), 2)


def s2tp(ra, dec, era=None, edec=None, raz=None, decz=None, unit="deg"):
    # Projection of spherical coordinates onto tangent plane:
    # "gnomonic" projection
    # Requires RA and Dec in degrees or unit set to "rad"
    # returns tangent plane coordinates in arcsec if unit=="deg", radians if
    #     unit =="rad"
    # tangent point selected automatically if not provided as raz, decz
    # advisable to let it do it itself unless one has a specific reason not to
    # errors will be propagated and returned if provided as era, edec
    # errors must be provided in milliarcsec

    # TODO: implement check on _denomj, within acceptable values (look at the starlink routine)
    # it doesnt appear to be necessary to correct when frame crosses midnight

    # calculates the centre point of the frame if not provided by the user
    raz = raz if raz else np.mean([np.amax(ra),np.amin(ra)])
    decz = decz if decz else np.mean([np.amax(dec),np.amin(dec)])

    if unit == "deg":
        # convert to radians
        ra = radians(ra)
        raz = radians(raz)
        dec = radians(dec)
        decz = radians(decz)
    elif unit != "rad":
        raise ValueError("unit: {0} not 'deg' or 'rad'".format(unit))

    # Trig functions
    _sdecz = sin(decz)
    _sdecj = sin(dec)
    _cdecz = cos(decz)
    _cdecj = cos(dec)
    _radifj = ra - raz
    _sradifj = sin(_radifj)
    _cradifj = cos(_radifj)

    # Reciprocal of star vector length to tangent plane
    _denomj = _sdecj * _sdecz + _cdecj * _cdecz * _cradifj

    # Compute tangent plane coordinates
    xi = _cdecj * _sradifj / _denomj
    xn = (_sdecj * _cdecz - _cdecj * _sdecz * _cradifj) / _denomj

    # convert to arcsec if necessary
    if unit == "deg":
        xi = np.degrees(xi)*3600.
        xn = np.degrees(xn)*3600.

    # a = ra of target
    # d = dec of target
    # m = ra of tp
    # n = dec of tp

    if era != None and edec != None:
        # xi derivatives
        dxi_da = (
                  (cos2(dec) * _cdecz * sin2(_radifj) /
                      np.power((_cdecj * _cdecz * _cradifj + _sdecj * _sdecz), 2)) +
                  (_cdecj * _cradifj / (_cdecj * _cdecz * _cradifj + _sdecj * _sdecz))
                 )
        dxi_dd = (
                  -(_sdecj * _sradifj / (_cdecj * _cdecz * _cradifj + _sdecj * _sdecz)) -
                  ((_cdecj * _sradifj * (_cdecj * _sdecz - _sdecj * _cdecz * _cradifj)) /
                      np.power((_cdecj * _cdecz * _cradifj + _sdecj * _sdecz), 2))
                 )

        # xn derivatives
        dxn_da = (
                  (_cdecj * _sdecz * _sradifj /
                      (_cdecj * _cdecz * _cradifj + _sdecj * _sdecz)) +
                  (_cdecj * _cdecz * _sradifj * (_sdecj * _cdecz - _cdecj * _sdecz * _cradifj) /
                      np.power(_cdecj * _cdecz * _cradifj + _sdecj * _sdecz, 2))
                 )
        dxn_dd = (
                  ((_sdecj * _sdecz * _cradifj + _cdecj * _cdecz) /
                      (_cdecj * _cdecz * _cradifj + _sdecj * _sdecz)) -
                  (((_cdecj * _sdecz - _sdecj * _cdecz * _cradifj) * (_sdecj * _cdecz - _cdecj * _sdecz * _cradifj)) /
                      np.power(_cdecj * _cdecz * _cradifj + _sdecj * _sdecz, 2))
                 )

        # standard error propagation
        exi = np.hypot(dxi_da * era, dxi_dd * edec)
        exn = np.hypot(dxn_da * era, dxn_dd * edec)
        
        if unit == "deg":
            raz, decz = np.degrees(raz), np.degrees(decz)    
        return xi, xn, exi, exn, (raz, decz)
    else:
        if unit == "deg":
            raz, decz = np.degrees(raz), np.degrees(decz)    
        return xi, xn, (raz, decz)

'''
def R(time):
    # R is the position vector of Earth
    # requires an astropy Time object
    # returns R, a three component array describing the position vector of Earth
    # in au
    # see eq. 12.25 in 'Spherical Astronomy' by R.M. Green

    # number of days since noon on 1/1/2000
    n = time - Time('2000-01-01T12:00:00', format='isot', scale='utc')

    # Obliquity of the ecliptic
    e = 23.439 - 0.0000004*n.jd
    e = radians(e)

    # mean longitude of the sun
    L = (280.460 + 0.9856474*n.jd) % 360.0
    # mean anomaly of the sun
    g = (357.528 + 0.9856003*n.jd) % 360.0
    # ecliptic longitude of the sun
    l = L + 1.915*sin(radians(g)) + 0.020*sin(radians(2*g))
    l = radians(l)

    return np.array([
                     -cos(l),
                     -sin(l) * cos(e),
                     -sin(l) * sin(e)
                    ]).T
'''


def W(a):
    # local west unit vector
    return np.array([
                     sin(a),
                     -cos(a),
                     0
                    ])


def N(a,d):
    # local north unit vector
    return np.array([
                     cos(a) * sin(d),
                     sin(a) * sin(d), # TODO: this is wrong according to green p339, should be negative
                     -cos(d)
                    ])


def T_pos(t, t0, coords, xi0, xn0, pmxi, pmxn, plx):
    # Returns a tangent plane position (xi, xn) at time 't'
    # Requirements:
    #   - astropy time objects 't' and 't0', the output times and start
    #     time respectively
    #   - astropy coordinate object describing the approximate position of the
    #     target
    #   - tangent point coordinates at time 't0'
    #   - proper motion of the target
    #   - parallax of the target
    # The units of xi0, xn0, pmxi, pmxn and plx are semi-arbitrary. Angular
    # units must simply all be the same. For proper motion (pmxi, pmxn) they
    # must be [angular] per year.

    # coords in radians
    a = coords.ra.rad
    d = coords.dec.rad

    # parallax factors in xi and xn
    Fxi = (R(t).dot(W(a)))
    Fxn = (R(t).dot(N(a,d)))

    # time difference in years
    dt = t.jyear - t0.jyear

    # coordinates at times in 't'
    xix = xi0 + dt * pmxi + plx * Fxi
    xnx = xn0 + dt * pmxn + plx * Fxn

    # format output
    if isinstance(xix, float):
        return np.array([xix, xnx])
    else:
        return np.array(zip(xix, xnx))


if __name__=="__main__":
    t0 = Time(2012.0, format='jyear')
    t = Time([2010.0, 2010.5, 2011.3, 2015.0], format='jyear')
    coord = SkyCoord(ra=180., dec=45., unit='deg', frame='icrs')
    print(T_pos(t, t0, coord, 0, 0, 2, 1, 0.5))
