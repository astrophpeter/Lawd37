########################
# A class for tracking the
# position of stars with
# propermotions
# @author P.McGill
#######################

from astropy.time import Time
import uncertainties.umath as np
from uncertainties import ufloat
from astropy.coordinates import SkyCoord

import parallax as p
import numpy as numpy
import microlens as m

class skyobj(object):
	
	MAS_TO_DEG = 1.0 / (3600.0 * 1000.0)
	DEG_TO_RAD = numpy.pi / 180.0

	def __init__(self,id,ra,dec,pmra,pmdec,epoch,parallax=None,Gmag=None):
		
		self._id = id
		self._ra_0 = ra
		self._dec_0 = dec
		self._pmra = pmra
		self._pmdec = pmdec
		self._epoch_0 = epoch
		self._parallax = parallax
		self._Gmag = Gmag


	def getRaDec(self,epoch):

		
		epoch_0 = Time(self._epoch_0,format='jyear',scale='tcb')
		epoch = Time(epoch,format='jyear',scale='tcb')
		
		dt = epoch.jyear - epoch_0.jyear

		if self._parallax == None:
			
			dec_0_rad = self._dec_0 * self.DEG_TO_RAD
			
			dec = self._dec_0 + dt * self._pmdec * self.MAS_TO_DEG
			ra = self._ra_0 + dt * (self._pmra / np.cos(dec_0_rad)) * self.MAS_TO_DEG
			
			return numpy.array([ra,dec])
		else:

			parallaxAU = p.parallax_to_au(self._parallax)

			R_0 = p.angular_to_cartesian(self._ra_0,self._dec_0) * parallaxAU

			_m = p.get_motion_vector(self._ra_0,self._dec_0,
							self._pmra,self._pmdec,self._parallax)
			
			R_earth = p.get_earth_observer_vector(epoch)

			R = R_0 + (_m* dt) - R_earth

			#This was used to test against astropy
			#Rcoord = SkyCoord((R[0]).n,(R[1]).n,(R[2]).n, frame='icrs', unit='AU', representation='cartesian')
			#Rcoord.representation = 'spherical'
			#return numpy.array([Rcoord.ra.degree,Rcoord.dec.degree])
			
			return p.cartesian_to_angular(R[0],R[1],R[2])



	def getSep(self,epoch,other):
		"""Get angular Separation with another skyobj,
		in mas
		"""

		pos1 = self.getRaDec(epoch)
		pos2 = other.getRaDec(epoch)

		ra1 = pos1[0]
		dec1 = pos1[1]
	
		ra2 = pos2[0]
		dec2 = pos2[1]
		

		#This was used to test against astropy
		#coord1 = SkyCoord(ra=ra1, dec=dec1, unit='deg', frame='icrs')
		#coord2 = SkyCoord(ra=ra2, dec=dec2, unit='deg', frame='icrs')
		#return coord1.separation(coord2).arcsec * 1000.0

		return m.get_angular_sep(ra1,dec1,ra2,dec2) 

					
	def getSepNum(self,epoch,other):
		"""Get angular Separation with another skyobj,
                in mas number not ufloat only for
		 minimization purposes
                """

		sep = self.getSep(epoch,other)

		
		if isinstance(sep,float):
			return sep
		else:
			return (sep).n
	
	def getPmAngle(self):

		angle = numpy.arctan2(self._pmdec,self._pmra)

		if angle  < 0:
			return (2*numpy.pi) + angle

		else: 
			return angle			
	
				
