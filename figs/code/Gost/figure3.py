
import numpy as np
from astropy.coordinates import SkyCoord
import vector_astrometry as va
import astrometry as a
from astropy.time import Time
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from astropy.io.fits import getdata
from matplotlib import gridspec
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
plt.style.use('idl.mplstyle')

# separation at a given julian year
def sepatt(jyear):
    time = Time(jyear, format='jyear')
    l = va.T_pos(time, t0, lcoord, lpmra, lpmdec, lplx, lvrad)
    s = va.T_pos(time, t0, scoord, spmra, spmdec, splx, svrad)
    return l.separation(s).arcsec * 1000.


# deflection calculation
def defl(sep):
    M = 0.61
    Dl = 4.63
    Oe = 90.25 * np.sqrt(M/Dl)
    u = sep/Oe
    return 0.5 * (np.sqrt(u**2 + 4) - u ) * Oe


# tangent plane coords of lens at given times
def tpcoordatT__l(time, raz, decz):
    l = va.T_pos(time, t0, lcoord, lpmra, lpmdec, lplx, lvrad)
    lxi, lxn, adz = a.s2tp(l.ra.deg, l.dec.deg, raz=raz, decz=decz, unit="deg")
    return lxi, lxn, adz
    
    
# tangent plane coords of source at given times
def tpcoordatT__s(time, raz, decz):
    s = va.T_pos(time, t0, scoord, spmra, spmdec, splx, svrad)
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
    
    # get min sep and time
    t_min = minimize_scalar(sepatt, method='bounded', bounds=[2010, 2030])
    print(t_min.x)
    print(sepatt(t_min.x))
    
    # times for trajectory plot
    tx = Time(np.arange(t_min.x-0.8, t_min.x+0.8, 1./365.2422), format='jyear')
    # times of Gaia transits
    
    t_gost_string = np.array([])

    for i in range(len(gost['TIME_UTC_NAME'])):
        print(str(gost['TIME_UTC_NAME'][i])[1:-2])
        t_gost_string = np.append(t_gost_string,str(gost['TIME_UTC_NAME'][i])[2:-1])


    #print(t_gost)
    t_gost = Time(t_gost_string, format='isot')
    
    # tangent plane coords at times for trajectory
    lxi, lxn, adz = tpcoordatT__l(tx, None, None)
    sxi, sxn, _ = tpcoordatT__s(tx, adz[0], adz[1])
    
    # tangent plane coords at Gaia transit times
    glxi, glxn, _ = tpcoordatT__l(t_gost, adz[0], adz[1])
    gsxi, gsxn, _ = tpcoordatT__s(t_gost, adz[0], adz[1])
    
    # deflections
    seps = sepatt(t_gost.jyear) # in mas
    defls = defl(seps) # in mas
    defl_angle = np.arctan2(gsxi-glxi, gsxn-glxn)
    
    # plot the deflections
    starty = 0.6
    defl_plot_factor = 0.1 
    # defl is in mas, plot axes are arcsec, multiplying defl by 0.1 means lines are 100x true length
    plt.text(0.0, 0.9, r"Deflection of Source ($\times$100), max %2.1f mas" % max(defls), horizontalalignment='center')
    for x, mag, angle in zip(glxi, defls, defl_angle):
        xoff = defl_plot_factor * mag * np.sin(angle)
        yoff = defl_plot_factor * mag * np.cos(angle)
        ax.arrow(x, starty, xoff, yoff, lw=1, length_includes_head=True, head_width=0.02, head_length=0.05, fc='k', ec='k')
        
    
    # plot the gaia scan directions
    starty = -0.6
    scan_plot_factor = 0.2
    plt.text(0.0, -0.9, r'Gaia scan direction', horizontalalignment='center')
    pos = 0
    for x, angle in zip(glxi, gost['scan_angle']):
        xoff = scan_plot_factor * np.sin(angle)
        yoff = scan_plot_factor * np.cos(angle)
        
        if Time(t_gost_string[pos],format='isot') < Time('2019-07-01',format='isot'):
            ax.arrow(x, starty, xoff, yoff, lw=1, length_includes_head=True, head_width=0.02, head_length=0.05, fc='k', ec='k')
        else:
            ax.arrow(x, starty, xoff, yoff, lw=1, length_includes_head=True, head_width=0.02, head_length=0.05, fc='grey', ec='grey')
        pos = pos +1

    plt.minorticks_on()

    plt.tick_params(axis='y', direction='in',which='both',right=True)
    plt.tick_params(axis='x',direction='in',which='both',top=True)
    
    # lens trajectory
    plt.plot(lxi, lxn, label='Lawd 37',linestyle='--',color='b')
    #plt.text(-0.9, 0.2, 'Lens', )
    
    # source trajectory
    plt.plot(sxi, sxn, label='Source',color='r')
    #plt.text(0.3, 0.45, 'Source')
   

 
    # lens pos during Gaia scans
    plt.scatter(glxi, glxn, marker='+',color='black')

    plt.plot(-1.53,-0.04,marker='+',color='black')
    plt.text(-1.4,-0.06,r'Gaia Observations')
    plt.plot([-1.62,-1.45],[-0.15,-0.15],linestyle='--',color='b')
    plt.text(-1.4,-0.16, r'LAWD 37')
    plt.plot([-1.62,-1.45],[-0.24,-0.24],linestyle='--',color='r')
    plt.text(-1.4,-0.25,r'Source')
    
    # figure related gubbins
    plt.xlim(-2.0, 2.0)
    plt.ylim(-1.0, 1.0)
    #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #       ncol=2, mode="expand", borderaxespad=0.)
    plt.xlabel(r'Relative $\alpha\cos\delta$ [arcseconds]')
    plt.ylabel(r'Relative $\delta$ [arcseconds]')
    plt.savefig('./fig3.eps',bbox_inches='tight',figsize=(8,12))
    plt.savefig('./fig3.png', dpi=600, bbox_inches='tight')
    plt.show()
    
    
    
    
    
    
    
    
