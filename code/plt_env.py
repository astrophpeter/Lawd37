from astropy.coordinates import SkyCoord
import numpy as np
import aplpy
from astropy.wcs import WCS
from astropy.io import fits, ascii
from astroquery.skyview import SkyView
from skyobj import skyobj
from astropy.time import Time

def datetime_to_jyTCB(date):
        """
        Returns the value of a time in Julian years from a date string
        of the format 'YYYY-MM-DD' scaled in (TCB)
        Args:
        date (string) : string format of the date 'YYYY-MM-DD'
        Return
        ulianYear (double) : the date in julian years

        Units Tests:
        #>>> math.isclose(unittestinglens.datetime_to_jyTCB('2009-01-01 00:00:00'),2009.0)
        True
        #>>> math.isclose(unittestinglens.datetime_to_jyTCB('2016-07-02'),2016.5)
        True
        """

        time = Time(date,scale='tcb')
        time.format = 'jyear'

        return time.value


def plt_lens_env(lens,source,event):
        """
        Creates and saves a plot of the source and lens
        environment. Usea DSS sky cut out images. In the
        plot the lens blue with trajectory in green. The
        source is identified in red. It is assumed source
        is stationary.

        Ags:
                lens (lens object): the lens star

                sourceRa (double): Source Right acession
                                   [Degrees]

                sourceDec (double) : Source Declination
                                    [Degrees]

        Returns:
                file (.png) : Output file saved with the id	
                      of the lens.
        """
        #define cutout image size [arsec]
        imsize = 2


        sourceCoords = source.getRaDec(2015.0)
        sourceRa = sourceCoords[0]
        sourceDec = sourceCoords[1]

        params = {'ra': sourceRa,'de': sourceDec,'imsize': imsize }

        #DSS Search
        url = "http://stdatu.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&r={ra}&d={de}&e=J2000&h={imsize}&w={imsize}&f=fits".format(**params)
        hdulist = fits.open(url)

        print(url)

        #Find the time the image was taken YYYY-MM-DD
        timeString  = hdulist[0].header['DATE-OBS'][:10]

        #Find the time the image was taken YYYY-MM-DD
        timeString  = hdulist[0].header['DATE-OBS'][:10]

        # Get the position of the lens at the time the image was taken
        time = datetime_to_jyTCB(timeString)
        coord = lens.getRaDec(time)

        #Get ra and dec of lens at the image time (raCos(dec) is converted into ra)
        raLensimag = coord[0]
        decLensimag = coord[1]

        # Find the postition of the lens in the future so its trajectory can be plotted
        timeend = datetime_to_jyTCB('2040-01-01')
        coordend = lens.getRaDec(timeend)
        raLensimagend = coordend[0]
        decLensimagend = coordend[1]

        # Find the postition of the lens at gaia epoch 2015.0
        coordGaiaEpoch = lens.getRaDec(2015.0)
        raLensimagGaiaEpoch = coordGaiaEpoch[0]
        decLensimagGaiaEpoch = coordGaiaEpoch[1]

        # Find the line of the lens trajectory between the time when the image was taken to
        # to some time in the future.
        line = np.array([[raLensimag,raLensimagend],[decLensimag,decLensimagend]])

        # Plot the image, with reverse gray color map.
        fig = aplpy.FITSFigure(hdulist[0])
        fig.show_colorscale(cmap='gray_r')

        #Plot the positions of the source, lens and lens tragectory on the image
        fig.show_markers(sourceRa,sourceDec,marker='o',edgecolor='r',label='source')
        fig.show_markers(raLensimag,decLensimag,marker='o',edgecolor='b',label='lens at image time (1989-11-22)')
        fig.show_lines([line],color='g',linestyle='--',label='lens-trajectory')

        #PLot gaia source positions also in the image
        #pos = get_gaia_source_pos(sourceRa,sourceDec,2)
        #fig.show_markers(pos[0],pos[1],marker='v',edgecolor='magenta',label='Gaia source')

        #Plot lens at gaia epoch on the image
        #fig.show_markers(raLensimagGaiaEpoch,decLensimagGaiaEpoch,marker='>',edgecolor='b')

        #Find Time and distance of Cloeset approach
        timeCl = event.get_time_of_minSep()
        distCl = event.get_min_sep()
        maxShift = event.get_max_resolved_centroid_shift()

        # get some info about the lens system to display on the plot
        Diststr = 'Distance of closest Approach: ' + str("%.2f" % distCl) + ' [mas]'
        Timestr = 'Time of closest Approach: ' + str("%.2f" %timeCl) + ' [Jyr]'
        Shiftstr = 'Max shift: ' + str("%.4f" %maxShift) + ' [mas]'

	#Lensstr = 'TGAS Lens Id (Blue): ' + str(lens.getId())
        #Sourcestr = 'PPMXL Source (Red) Id:' + str(sourceId)
        #SourceMag = 'PPMXL source mag [j,h,k,b1,b2,r1,r2]: '
        #LensMag = 'TGAS Lens Gmag: ' + str("%.2f" %l

        # Add some extra info about the source and lens to the plot.
        fig.add_label(0.75,0.98,Diststr,relative=True)
        fig.add_label(0.75,0.94,Timestr,relative=True)
        fig.add_label(0.75,0.90,Shiftstr,relative=True)
	#fig.add_label(0.75,0.90,Lensstr,relative=True)
        #fig.add_label(0.75,0.86,Sourcestr,relative=True)
        #fig.add_label(0.80,0.82,SourceMag,relative=True)
        #fig.add_label(0.70,0.78,str(sourceMag),relative=True)
        #fig.add_label(0.85,0.74,LensMag,relative=True)


        filename = 'candSoft/TGAS_' + str(lens._id) + '_' + str(source._id) + '.png'
        fig.save(filename,dpi=200)
