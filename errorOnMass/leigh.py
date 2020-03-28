#!/usr/bin/env python

# define a few functions that I use regularly

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time


def crappyhist(a, bins):
	'''Draws a crappy text-mode histogram of an array'''
	import string
	from math import log10
	
	h,b = np.histogram(a, bins)
	
	for i in range (0, bins-1):
	#	print string.rjust(`b[i]`, 7)[:int(log10(
#				   np.amax(b)))+5], '| ', '#'*int(70*h[i-1]/np.amax(h))
	
	###print string.rjust(`b[bins]`, 7)[:int(log10(np.amax(b)))+5] 
		Pass



def eq2gal(a, d, mua, mud, emua=None, emud=None, errs=False):
	coords = apsc(a, d)
	#https://arxiv.org/pdf/1306.2945.pdf
	a, d = np.radians(a), np.radians(d)
	ag = np.radians(192.85948)
	dg = np.radians(27.12825)
	lo = np.radians(32.93192)

	C1 = np.sin(dg)*np.cos(d) - np.cos(dg)*np.sin(d)*np.cos(a-ag)
	C2 = np.cos(dg)*np.sin(a-ag)

	f = 1.0/np.hypot(C1,C2)
	mulb = f * np.array([C1*mua + C2*mud, -C2*mua + C1*mud])
	if errs:
		emulb = f * np.array([C1*emua + C2*emud, -C2*emua + C1*emud])
		emul = np.hypot( (f*C1*emua), (f*C2*emud) )
		emub = np.hypot( (f*-C2*emua), (f*C1*emud) )
		return coords.galactic.l.deg, coords.galactic.b.deg, mulb[0], mulb[1], emul, emub

	return coords.galactic.l.deg, coords.galactic.b.deg, mulb[0], mulb[1]



def apsc(ra, dec):
	# return a position as astropy sky coordinate
	# requires:
	#	RA and Dec in decimal degrees (these can be 1d arrays)
	return SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')


def apt(t, fmat=None):
	# return a time as an astropy time
	# requires:
	#	time
	#	format
	if not fmat:
		print("No format given for time, must be one of:")
		print(list(Time.FORMATS))
		return
	return Time(t, format=fmat)


def pmpos(coord0, t0, pmra, pmde, t):
	# return a position at a given epoch
	# taking pm into account
	# requires:
	#	astropy Sky Coord instance for initial epoch
	#	astropy time instance for initial epoch
	#	proper motion in RA and Dec (mas)
	#	astropy time instance for final epoch
	tdiff = (t.decimalyear-t0.decimalyear)
	ra = coord0.ra.deg + (pmra/3.6E6)/np.cos(np.radians(coord0.dec.deg))*tdiff
	de = coord0.dec.deg + (pmde/3.6E6)*tdiff
	return apsc(ra, de)

def absMag(plx, eplx, mag, emag):
	# return an absolute magnitude and error
	# requires:
	#	parallax and error (in mas)
	#	magnitude and error
	Mag = mag + 5.0*(1.0+np.log10(plx/1000.0))
	eMag = np.sqrt(
		emag**2 +
		((5.0*(eplx/1000.0))/((plx/1000.0)*np.log(10.0)))**2)
	return Mag, eMag



def plxDist(plx, eplx):
	# return a distance and error
	# requires:
	#	parallax and error (in mas)
	Dist = 1000.0/plx
	eDist = np.sqrt(((-1000.0*eplx)/(plx**2))**2)
	return Dist, eDist



def doMedian(data):
	# return a median and standard error
	# requires:
	#	array of data
	data = np.float_(data)
	avg = np.median(data)
	eavg = 1.253 * np.std(data)/np.sqrt(len(data))
	return avg, eavg



def doMean(data):
	# return a mean and standard error
	# requires:
	#	array of data
	data = np.float_(data)
	avg = np.mean(data)
	eavg = np.std(data)/np.sqrt(len(data))
	return avg, eavg



def projSep(AngSep, eAngSep, Plx, ePlx):
	# returns a projected separation in AU
	# requires:
	#	angular separation and error in arcsec
	#	parallax and error in mas
	Plx /= 1000.0
	ePlx /= 1000.0
	ProjSep = AngSep/Plx
	eProjSep = np.hypot(
		(eAngSep/Plx),
		((-AngSep*ePlx)/(Plx**2)))
	return ProjSep, eProjSep



def orbVel(M, a, r):
	# returns orbital velocity in km/s
	# requires:
	# 	system mass (solar masses)
	# 	orbital semi-major axis (in AU)
	# 	orbital radius (in AU)
	G = 6.67384E-11
	Msun = 1.9891E+30
	AU = 1.49597871E+11
	mu = G*M*Msun
	r, a = r*AU, a*AU
	V = np.sqrt(mu*((2/r)-(1/a)))
	return V/1000.0




def weightedMean(values, errors):
	# returns a weighted mean and associated uncertainty
	# requires:
	# 	an array of values
	# 	an array of equivalent errors (i.e. same length, same order)
	values = np.array(values)
	errors = np.array(errors)
	weights = 1.0/(errors**2)
	Avg = np.sum(values*weights)/np.sum(weights)
	eAvg = 1.0/(np.sqrt(np.sum(weights)))
	return Avg, eAvg
