from skyobj import skyobj
import microlens as m
from event import event
import numpy as numpy
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
plt.style.use('idl.mplstyle')
import uncertainties as np
from uncertainties import ufloat
from lawd37event import *
from astropy.time import Time

signal = numpy.array([])
distance = numpy.array([])
distance_err = numpy.array([])
signal_err = numpy.array([])
us = numpy.array([])
us_err = numpy.array([])

minTime = MicroEvent.get_time_of_minSep()
maxSignal = MicroEvent.get_max_resolved_centroid_shift()
ER = MicroEvent.get_einstein_R()
minDist = MicroEvent.get_min_sep()
minU = minDist / ER

print(minTime)
print(minDist)

times = numpy.linspace(minTime -0.5,minTime+0.5,num=250)

for t in times:
	signal = numpy.append(signal,(MicroEvent.get_resolved_centroid_shift_at_epoch(t)).n)
	u = MicroEvent.get_sep(t) / ER
	us = numpy.append(us,(u).n)
	
	distance = numpy.append(distance,(MicroEvent.get_sep(t)).n)
	distance_err = numpy.append(distance_err,(MicroEvent.get_sep(t)).s)
	
	us_err = numpy.append(us_err,(u).s)
	signal_err = numpy.append(signal_err,(MicroEvent.get_resolved_centroid_shift_at_epoch(t)).s)


#distance = distance / MicroEvent.get_einstein_R()


times = times - minTime

fig = plt.figure()
#set height ratios for sublots
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2])

#the fisrt subplot
ax0 = plt.subplot(gs[0])
ax0.plot(times,signal,color='black',linewidth=0.3)
ax0.fill_between(times,signal-signal_err,signal+signal_err,facecolor='gray',interpolate=True,alpha=0.5)
ax0.fill_between(times,signal-signal_err,signal-2*signal_err,facecolor='gray',interpolate=True,alpha=0.25)
ax0.fill_between(times,signal+signal_err,signal+2*signal_err,facecolor='gray',interpolate=True,alpha=0.25)
ax0.scatter(0.0,(maxSignal).n,marker='o',color='red',alpha=1.0)

# + '\pm' + str(round((maxSignal).s,1))+'$')


#the second subplot
# shared axis X
ax1 = plt.subplot(gs[1], sharex = ax0)
ax1.plot(times,distance, color='black',linewidth=0.3)
ax1.fill_between(times,distance-distance_err,distance+distance_err,facecolor='gray',interpolate=True,alpha=0.5)
ax1.fill_between(times,distance-distance_err,distance-2*distance_err,facecolor='gray',interpolate=True,alpha=0.25)
ax1.fill_between(times,distance+distance_err,distance+2*distance_err,facecolor='gray',interpolate=True,alpha=0.25)
ax1.scatter(0.0,(minDist).n,marker='s',color='b',alpha=1.0)

ax12 = ax1.twinx()
ax12.plot(times,us,color='white',linewidth=0.0)
ax12.set_ylabel(r'$u$')

#ax13 = ax0.twiny()
#ax13.set_xlim(times.min()-0.05,times.max()+0.05)

#Getting human readable times
humanTime = Time([-0.4,0.2,0.0,0.2,0.4] + minTime,format='jyear')
humanTime = humanTime.isot

#ax13.set_xticks((-0.5,0.0,0.5),
#           (str(humanTime[0]), str(humanTime[125]),str(humanTime[249])))




#ax0.xaxis.set_major_formatter(plt.NullFormatter())

ax1.set_ylabel(r'$\Delta\theta$ [mas]')
ax0.set_ylabel(r'$\delta\theta$ [mas]')
ax1.set_xlabel('Time [yyyy-mm-dd]')
ax1.set_ylim(0,1000)
ax1.axhline(y=100,linestyle='--',color='r')

#yticks = ax1.yaxis.get_major_ticks()
#yticks[-1].label1.set_visible(False)

ax1.minorticks_on()
ax0.minorticks_on()
ax12.minorticks_on()
#ax13.minorticks_on()
ax1.tick_params(axis='y', direction='in',which='both',top=True)
ax1.tick_params(axis='x',direction='in',which='both',top=True)
#ax13.tick_params(axis='x',direction='in',which='both',top=True)
ax12.tick_params(axis='y',direction='in',which='both')
ax12.set_ylim(0,30)

ax0.tick_params(axis='x', direction='in',which='both',right=True,top=True)
ax0.tick_params(axis='y',direction='in',which='both',right=True)

# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)


#Plot Legends
size = 'small'
dy = 0.13
x0 = 0.1
y0 = 3.0
dx = 0.025


ax0.text(x0+dx,y0,r'$\delta\theta_{max}=' +str(round((maxSignal).n,1))+'\pm'+str(round((maxSignal).s,1))+'$ mas',size=size)
ax0.scatter(x0,y0+0.5*dy,marker='o',color='r')

ax0.text(x0+dx,y0-2*dy,r'$u_{min}=' + str(round((minU).n,1)) + '\pm'+ str(round((minU).s,1))+ '$',size=size)
ax0.scatter(x0,y0-1.5*dy,marker='s',color='b')

ax0.text(x0+dx,y0-dy,r'$\Delta\theta_{min}=' + str(int(round((minDist).n,-1))) + '\pm' + str(int(round((minDist).s,-1))) + '$ mas',size=size)
ax0.scatter(x0,y0-0.5*dy,marker='s',color='b')

ax0.plot([-0.52,-0.46],[y0,y0],linestyle='--',color='r')
ax0.text(-0.45,y0-dy,'Gaia\'s source resolution \n limit $\sim$100 mas',size=size)
gs.update()
fig.set_figwidth(5)
fig.set_figheight(6.5)

yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)

#ax0.set_xticks([-0.5,0.0,0.5],
#           [str(humanTime[0]), str(humanTime[125]),str(humanTime[249])])

#ax0.xaxis.set_ticks_position('top')
ax0.xaxis.set_label_position('top')
#plt.tight_layout()

fig.canvas.draw()

labels = [item.get_text() for item in ax1.get_xticklabels()]
labels[1] = humanTime[0][:10]
labels[3] = humanTime[2][:10]
labels[5] = humanTime[4][:10]

ax1.set_xticklabels(labels)


plt.savefig('test.pdf',figsize=(10,12))


