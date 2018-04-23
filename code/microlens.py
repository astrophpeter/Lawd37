import uncertainties.umath as np
import numpy as numpy


deg2rad = numpy.pi / 180.0

def get_angular_sep1(ra1,dec1,ra2,dec2):

	"""Calculates the angular separation of two points on the
        sphere using great circles distance formula.

        returns angular separation in [mas]
        """

	ra1 = ra1 * deg2rad
	ra2 = ra2 * deg2rad

	dec1 = dec1 * deg2rad
	dec2 = dec2 * deg2rad

	dist = np.acos(np.sin(dec1)*np.sin(dec2) + (np.cos(dec1)*np.cos(dec2) * np.cos(ra1-ra2)))
	
	
	return (dist / deg2rad) * 3600.0 * 1000.0

def get_angular_sep(ra1,dec1,ra2,dec2):
	"""
	implements the vincinerty distance formular in the case
	of a sphere.
	"""

	ra1 = ra1 * deg2rad
	ra2 = ra2 * deg2rad


	dec1 = dec1 * deg2rad
	dec2 = dec2 * deg2rad

	top = np.sqrt((np.cos(dec2)*np.sin(ra2-ra1))**(2)+(np.cos(dec1)*np.sin(dec2)-np.sin(dec1)*np.cos(dec2)*np.cos(ra2-ra1))**(2))

	bottom = (np.sin(dec1)*np.sin(dec2))+(np.cos(dec1)*np.cos(dec2)*np.cos(ra2-ra1))

	dist = np.atan2(top,bottom)

	return (dist / deg2rad) * 3600.0 * 1000.0

def get_einstein_R(lensMass,lensDist,sourceDist=None):
	"""Calculates the einstein radius of a lens and
	source system.


	Args:
	   lensMass (float) : Mass of the foreground lens object [Msol]

           lensMass Error (float) :  Error im the Mass of the foreground lens
                object [Msol]

           lensDist (float) : Distance to the lens [pc]

           lensDist Error (float) : Error in the Distance to the lens [pc]

           SourceDist (float,optional) : Distance to the source.
                                         defaults to None. If none,
                                         caculation will assume sourceDist
                                         to be at infinity. [pc]


        Returns:
           einsteinR (float) : The Einstein radius [mas]

	"""


	if sourceDist is None:
		return (90.2 * np.sqrt(lensMass / lensDist))

	else:
		return (90.2 * np.sqrt((lensMass / lensDist) * (1- lensDist / sourceDist)))


def get_centroid_shift_unresolved(lensMass,lensDist,Sep,sourceDist=None,
                                lensMag=None,sourceMag=None):
	"""Calculates the expected astrometric centriod shift
	of a lens source system.

	Args:
	   lensMass (float) : lensMass (float) : Mass of the foreground
                              lens object [Msol]

           lensDist (float) : Distance to the lens
                              [pc]

           SourceDist (float,optional) : Distance to the source.
                                         defaults to None. If none,
                                         caculation will assume sourceDist
                                         to be at infinity. [pc]

           Sep (float) : Angular separation between the lens
                            abd source [mas]

           lensMag (float,optional) : Magnitude of the lens. Defaults to None.
                                      If none will assume the dark lens equation.

           sourceMag (float,optional) : Magnitude of the source. Defaults to None.
                                        If none will assume the dark lens eqution.

         Returns:
            centriodShift (float) : The expected centriod shift [mas]
	"""

	EinsteinR = get_einstein_R(lensMass,lensDist,sourceDist=sourceDist)
	mu = Sep / EinsteinR

	if lensMag is None and sourceMag is None:

		return (mu * EinsteinR) / (mu**2 +2)

	else:
		lumFactor = 1+(10** ((lensMag - sourceMag) / (-2.5)))
		return (mu * EinsteinR) / ((mu**2 + 2)*lumFactor)


def get_centriod_shift_resolved(lensMass,lensDist,sep,sourceDist=None):
	"""Calculates the centroid shift of the source position
	due to the brightest image.

	See Eq (2) of Sahu et al (2017)

	Args:
           lensMass (float) : lensMass (float) : Mass of the foreground
                              lens object [Msol]

           lensDist (float) : Distance to the lens
                              [pc]

           SourceDist (float,optional) : Distance to the source.
                                         defaults to None. If none,
                                         caculation will assume sourceDist
                                         to be at infinity. [pc]
	"""

	einsteinR = get_einstein_R(lensMass,lensDist,sourceDist=sourceDist)
	u = sep / einsteinR

	return get_major_img(lensMass,lensDist,sep,sourceDist=sourceDist) - (u*einsteinR)


def get_major_img(lensMass,lensDist,sep,sourceDist=None):

	einsteinR = get_einstein_R(lensMass,lensDist,sourceDist=sourceDist)
	u = sep / einsteinR

	return 0.5 * (np.sqrt(u**2 +4)+u) * einsteinR


def get_mass_onemeasure(lensDist,sep,shift):

	return ((1/90.2)**2)*(lensDist)*((shift))*((shift+sep))

def get_mass_fractional_err(lensDist,EnstienR,parallax,astrometric_err):

	return (astrometric_err /1000.0)*np.sqrt(4*(lensDist*1000.0/(EnstienR))**2+(1000.0/parallax)**2)


