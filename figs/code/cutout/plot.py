import aplpy
import matplotlib.pyplot as mpl
import sys
sys.path.append('../../skyobj/')
import numpy as np
from skyobj import skyobj
from astropy.io import fits,ascii
from astropy.time import Time
import pyfits

from matplotlib.path import Path 
verts = [(0.0, -1.0),(0.0, -0.3),(0.0, +1.0),(0.0, +0.3),(-1.0, 0.0),(-0.3, 0.0),(+1.0, 0.0),(+0.3, 0.0)]
codes = [Path.MOVETO,Path.LINETO,Path.MOVETO,Path.LINETO,Path.MOVETO,Path.LINETO,Path.MOVETO,Path.LINETO,]
path = Path(verts, codes)


fig = mpl.figure(figsize=(11,20))

lens = skyobj(id=1,ra=176.454907296219, dec=-64.842957135494, pmra=2662.03572627, pmdec=-345.18255501, parallax=215.78,epoch=2015.0)
source = skyobj(id=2,ra=176.46360456073, dec=-64.8432977866831, pmra=-14, pmdec=-2,epoch=2015.0)


times = np.linspace(2012,2030,num=150)
ras = np.array([])
decs = np.array([])
for t in times:
	ra,dec = lens.getRaDec(t)
	ras = np.append(ras,ra)
	decs = np.append(decs,dec)

#print(ras)
#print(decs)

#decaphdu = fits.open('decaps(1).fits')

# Launch APLpy figure of image
imgL = aplpy.FITSFigure('dssred.fits',figure=fig,subplot=(2,1,1))
imgR = aplpy.FITSFigure('noartifactDecaps.fits',figure=fig,subplot=(2,1,2))
imgInset = aplpy.FITSFigure('decapSupersZoomTight.fits',figure=fig,subplot=[0.1,0.38,0.52,0.08])
#imgInset = aplpy.FITSFigure('decapSupersZoomTight.fits',figure=fig,subplot=[0.52,0.73,0.33,0.17])
# Apply grayscale mapping of image


hdulist = fits.open('dssred.fits')

timestring = hdulist[0].header['DATE-OBS'][:10]
time = Time(timestring,scale='tcb')
time.format = 'jyear'
timedecimal = time.value 

print(timedecimal)

# Or apply a different stretch to the image

imgL.show_colorscale(cmap='gray_r',stretch='log')
imgL.add_label(0.9,0.96,'Epoch {:.0f}'.format(timedecimal),relative=True,size='x-large')

#get lens position at time of image and future
lensRaIm,lensDecIm = lens.getRaDecNoPlx(timedecimal)
lensRaF,lensDecF = lens.getRaDecNoPlx(2015.0)
lensRaFF,lensDecFF = lens.getRaDecNoPlx(2030.0)

lensLine = np.array([[lensRaIm,lensRaF],[lensDecIm,lensDecF]])

sourceRaIm,sourceDecIm = source.getRaDecNoPlx(timedecimal)
sourceRaF,sourceDecF = source.getRaDecNoPlx(2040.0)
sourceLine = np.array([[sourceRaIm,sourceRaF],[sourceDecIm,sourceDecF]])

#get time of closest approach
minTime = lens.getMinTime(source)
lensRaC,lensDecC = lens.getRaDec(minTime)

#get source position at time of image and future
sourceRaC,sourceDecC = source.getRaDecNoPlx(minTime)

imgL.show_lines([lensLine],color='b',linestyle='--')
imgL.show_markers(lensRaIm,lensDecIm,edgecolor='b',marker='o',s=5000.0)
imgL.show_arrows(lensRaIm,lensDecIm,lensRaF-lensRaIm,lensDecF-lensDecIm,edgecolor='b',width=0.5,head_width=2,head_length=2)
imgL.add_label(lensRaIm,lensDecIm - 0.003,'Lawd 37',size='x-large',color='b')
imgL.add_label(sourceRaC,sourceDecC - 0.003,'Gaia Source \n 5332606346467258496',color='r',size='large')
#imgL.show_markers(lensRaIm,lensDecIm,layer='markers',edgecolor='b',facecolor='none',marker=path,s=10000.0,alpha=0.5,linewidths=2)
#imgL.show_markers(lensRaC,lensDecC,edgecolor='b',marker='*')

imgL.show_markers(sourceRaIm,sourceDecIm,color='r',marker='o',s=500.0)
imgL.ticks.hide()
imgL.frame.set_linewidth(0)

lensRaImD,lensDecImD = lens.getRaDec(2016.5)
sourceRaImD,sourceDecImD = source.getRaDecNoPlx(2016.3)

#imgR.show_lines([lensLine],color='b',linestyle='--')
imgR.show_markers(lensRaImD,lensDecImD,edgecolor='b',marker='o',s=2000.0)


#imgR.show_markers(lensRaC,lensDecC,edgecolor='b',marker='*')

imgR.show_markers(sourceRaC,sourceDecC,color='r',marker='o',s=200.0)
imgR.show_arrows(lensRaImD,lensDecImD,lensRaFF-lensRaImD,lensDecFF-lensDecImD,edgecolor='b',width=0.5,head_width=7,head_length=7)
# Or apply a different stretch to the image
imgR.show_colorscale(cmap='gray_r',stretch='arcsinh')
imgR.add_label(0.9,0.96,'Epoch 2016',relative=True,size='x-large')
imgR.axis_labels.hide()
imgR.ticks.hide()
imgR.frame.set_linewidth(0)

#ImgR.show_regions('ds9.reg')

imgInset.show_colorscale(cmap='gray_r',stretch='arcsinh')
#imgInset.show_markers(lensRaImD,lensDecImD,edgecolor='b',marker='o',s=100)

imgInset.show_markers(lensRaC,lensDecC,edgecolor='b',facecolor='b',marker='*',s=100)
imgInset.show_markers(sourceRaC,sourceDecC,edgecolor='r',facecolor='r',marker='*',s=100)
#imgInset.show_markers(sourceRaImD,sourceDecImD,color='r',marker='o')

imgInset.show_markers([ras],[decs],edgecolor='b',facecolor='b',marker='o',s=7.0)
#imgInset.show_lines([sourceLine],color='r',linestyle='--')

imgInset.tick_labels.hide()
imgInset.axis_labels.hide()
imgInset.ticks.hide()

imgL.tick_labels.hide()
imgL.axis_labels.hide()

imgR.ticks.set_color('black')
imgR.tick_labels.hide()
imgR.axis_labels.hide()
#ImgR.axis('off')


imgL.add_scalebar(0.013888888, "50 ''", color='black', corner='bottom right')
imgR.add_scalebar(0.013888888, "50 ''", color='black', corner='bottom right')

fig.tight_layout()
fig.savefig('ic348_basic.eps')

