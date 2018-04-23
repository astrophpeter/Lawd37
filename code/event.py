#####################
# A class to simulate 
# a microlensing event
# @author P.McGill
#####################


import numpy as numpy
import uncertainties.umath as np
from scipy.optimize import minimize_scalar
from skyobj import skyobj
from uncertainties import ufloat
import microlens as m

class event(object):

	def __init__(self,lens,source,lensMass,lensDist,sourceDist=None):

		self._lens = lens
		self._source = source
		self._lensMass = lensMass
		self._lensDist = lensDist
		self._sourceDist = sourceDist

	

	def get_time_of_minSep(self):
		"""Gets time of closest approach in Julian Years
		"""

		minimum = minimize_scalar(self._lens.getSepNum,
			args=(self._source),bounds=(2015,2025),method='bounded')
		return minimum.x

	def get_min_sep(self):
		"""Returns minimum Separation
		"""

		minTime = self.get_time_of_minSep()
		return self._lens.getSep(minTime,self._source)
	
	def get_sep(self,epoch):
		
		return self._lens.getSep(epoch,self._source)


	def get_einstein_R(self):
		"""Returns the Enstien Radius [mas]

		"""

		return m.get_einstein_R(self._lensMass,self._lensDist,sourceDist=self._sourceDist)



	def get_max_resolved_centroid_shift(self):
		"""Returns the centriod shift due to mmajor image only.

		"""

		MinSep = self.get_min_sep()		

		return m.get_centriod_shift_resolved(self._lensMass,
					self._lensDist,MinSep,sourceDist=self._sourceDist)



	def get_resolved_centroid_shift_at_epoch(self,epoch):
		"""Returns unresolved centriod shift at epoch [mas]
		"""

		sep = self._lens.getSep(epoch,self._source)
		
		return m.get_centriod_shift_resolved(self._lensMass,
                                        self._lensDist,sep,sourceDist=self._sourceDist)			


	def get_unresolved_centroid_shift_at_epoch(self,epoch):
		"""Returns mag of unresolved centriod shift at epoch

		"""
		sep = self.get_sep(epoch)
		ER = self.get_einstein_R()
		u = sep / ER
		return (ER *u)/(u**2+2)
		

	def get_max_unresolved_centroid_shift(self):

		MinSep = self.get_min_sep()
		return m.get_centroid_shift_unresolved(self._lensMass,self._lensDist,MinSep)


	def get_source_apparent_pos(self,epoch,enlargeFactor=1.0):

		lensPos = self._lens.getRaDec(epoch)
		sourcePos = self._source.getRaDec(epoch)
		
		shift_mag = self.get_unresolved_centroid_shift_at_epoch(epoch)
		shift_angle = np.atan2((sourcePos[0]-lensPos[0]) * np.cos(numpy.deg2rad(lensPos[1])),sourcePos[1]-lensPos[1])


		appSourceRa = sourcePos[0] + ((enlargeFactor * shift_mag * np.cos(shift_angle)) / (3600 * 1000))
		appSourceDec = sourcePos[1] + ((enlargeFactor * shift_mag * np.sin(shift_angle)) / (36000 * 1000))

		return numpy.array([appSourceRa,appSourceDec])
	
		
	def get_source_chi_squared_ra(self,epochStart,epochEnd,numSamps=10):

		chiS = 0;

		for i in numpy.linspace(epochStart,epochEnd,num=numSamps):
			appPos = self.get_source_apparent_pos(i)
			truePos = self._source.getRaDec(i)
			
			sep = ((appPos[1] - truePos[1]))**2
			#sep = ((appPos[0] - truePos[0])*numpy.cos(numpy.deg2rad(truePos[1])))**2
			chiS = chiS + (sep / abs(truePos[1]))  

		return chiS	
	
	def get_motion_vec_angle(self,epoch):
		"""return angle relative to constance lat"""


		SourcePos = self.get_source_apparent_pos(epoch)
		SourcePosdt = self.get_source_apparent_pos(epoch + 0.01)

		angle = numpy.arctan2(SourcePos[1]-SourcePosdt[1],numpy.cos(numpy.deg2rad(SourcePos[1]))*(SourcePos[0]-SourcePosdt[0]))

		if angle < 0:	
			return (2 * numpy.pi) + angle
		else:
			return angle


	def get_rough_pm(self,epoch_start,epoch_end):

		"""
		get rough proper motion between two points for the source due to lensing

		"""
			
		dt = epoch_end - epoch_start

		start = self.get_source_apparent_pos(epoch_start)
		end = self.get_source_apparent_pos(epoch_end)

		return numpy.array([(end[0] - start[0]) / dt, (end[1] - start[1]) / dt]) * 3600.0 * 1000.0 
