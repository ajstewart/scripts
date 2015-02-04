#!/usr/bin/python
#
#
# Author: Antonia Rowlinson
# E-mail: b.a.rowlinson@uva.nl
#
import sys
import tools
import tkp.utility.coordinates as coords
import numpy as np
import math
import os
import pytz
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
pylab.rcParams['legend.loc'] = 'best'
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties

###################### INITIAL SETUP STEPS ######################

if len(sys.argv) != 10:
    print 'python TraP_QC_diagnostics.py <database> <username> <password> <host> <port> <databaseType> <dataset_id> <skymodel> <telescope>'
    exit()
database = sys.argv[1]
username = sys.argv[2]
password = sys.argv[3]
host = sys.argv[4]
port = sys.argv[5]
databaseType = sys.argv[6]
dataset_id = str(sys.argv[7])
skymodel = str(sys.argv[8])
telescope = sys.argv[9]

if telescope != 'LOFAR' and telescope != 'MWA':
    print 'Elevation and sidereal time plots currently only support LOFAR and MWA.'
    print 'Edit line 31 of TraP_QC_flux.py to include your telescope name,'
    print 'and include the IRTF position of your telescope in the sidetime function'
    print 'given in tools.py.'
    print 'The elevation and sidereal time plots will otherwise be blank.'

###################### MAIN SCRIPT ######################

if not os.path.exists('ds'+str(dataset_id)+'_plotdata.txt'):

    # Extracting data from the TraP database into a text file
    if not os.path.exists('ds_'+dataset_id+'_sources_fluxQC.csv'):
        tools.get_data(database, username, password, host, port, databaseType, dataset_id)

    avg_flxrat=[]
    sources = tools.extr_src_data(dataset_id)
    sky_data={}
    flux_data=[]

    unique_pointings= np.unique([str(a[3])+','+str(a[2]) for a in sources])
    unique_frequencies = np.unique([int(float(a[8])/1e6) for a in sources])
    print unique_frequencies
    
    unique_sources = np.unique([a[10] for a in sources])
    
    skymodel_data={}
    
    if skymodel == 'auto':
        for pointing in unique_pointings:
            ptRA=pointing.split(',')[0]
            ptDec=pointing.split(',')[1]
            for frequency in unique_frequencies:
                if skymodel == 'auto':
                    skymodel_name=str(int(round(float(ptRA),0)))+'_'+str(int(round(float(ptDec),0)))+'.sky'
                    skymodel_key=str(int(round(float(ptRA),0)))+str(int(round(float(ptDec),0)))
                    search_radius=20.
                    flux_limit=0.0
                    if not os.path.isfile(skymodel_name):
                        os.system('gsm.py '+skymodel_name+' '+str(ptRA)+' '+str(ptDec)+' '+str(search_radius)+' '+str(flux_limit))
                    if skymodel_key not in sky_data.keys():
                        vlss_sources=[]
                        vlss_data=open(skymodel_name, 'r')
                        lines = iter(vlss_data)
                        line1=vlss_data.readline()
                        frq=float(line1.split("'")[1])/1e6
                        lines = iter(vlss_data)
                        lines.next()
                        lines.next()
                        for line in lines:
                            vlss_sources.append(line)
                        vlss_data.close()
                        sky_data[skymodel_key]=vlss_sources
                    skymodel_data[pointing+','+str(frequency)]=tools.extract_sky(sky_data[skymodel_key],float(frequency),frq)
    elif skymodel != 'auto' and skymodel != 'none':
        skymodel_data[skymodel]= tools.extr_skymodel(skymodel)

    source_data=[]
    fluxDPTs={}
    src_assoc={}
    avgs={}

    for x in sources:
        if str(x[11]) not in fluxDPTs.keys():
            fluxDPTs[str(x[11])]=[[float(x[6]),int(float(x[8])/1e6)]]
            if skymodel == 'auto':
                skymodel_name = str(x[3])+','+str(x[2])+','+str(int(float(x[8])/1e6))
                skymodel_flux=tools.source_assoc(skymodel_data[skymodel_name],x,float(x[9]))
            elif skymodel != 'auto' and skymodel != 'none':
                skymodel_flux=tools.source_assoc(skymodel_data[skymodel],x,float(x[9]))
            else:
                skymodel_flux=-1.
            src_assoc[str(x[11])]=skymodel_flux
        else:
            fluxDPTs[str(x[11])].append([float(x[6]),int(float(x[8])/1e6)])

    for freq in unique_frequencies:
        for srckey in fluxDPTs.keys():
            avgs[srckey+','+str(freq)]=np.average([x[0] for x in fluxDPTs[srckey] if x[1]==freq])

    if telescope == 'LOFAR' or telescope == 'MWA':
        source_data = [[float(x[16]),float(x[15]),float(x[6]),float(x[7]),avgs[str(x[11])+','+str(int(float(x[8])/1e6))],
                    coords.julian_date(datetime.strptime(x[13].split('.')[0], '%Y-%m-%d %H:%M:%S').replace(tzinfo=pytz.utc), modified=True),
                    int(float(x[8])/1e6),coords.angsep(float(x[16]),float(x[15]),float(x[3]),float(x[2]))/3600.,
                    tools.sidetime(x[16], x[15], coords.julian_date(datetime.strptime(x[13].split('.')[0], '%Y-%m-%d %H:%M:%S').replace(tzinfo=pytz.utc), modified=True), telescope)[1],
                    tools.sidetime(x[16], x[15], coords.julian_date(datetime.strptime(x[13].split('.')[0], '%Y-%m-%d %H:%M:%S').replace(tzinfo=pytz.utc), modified=True), telescope)[0],
                    src_assoc[str(x[11])], float(x[10]), x[14]] for x in sources]
    else:
        source_data = [[float(x[16]),float(x[15]),float(x[6]),float(x[7]),avgs[str(x[11])+','+str(int(float(x[8])/1e6))],
                    coords.julian_date(datetime.strptime(x[13].split('.')[0], '%Y-%m-%d %H:%M:%S').replace(tzinfo=pytz.utc), modified=True),
                    int(float(x[8])/1e6),coords.angsep(float(x[16]),float(x[15]),float(x[3]),float(x[2]))/3600.,
                    -1,
                    -1,
                    src_assoc[str(x[11])], float(x[10]), x[14]] for x in sources]
    
    output3 = open('ds'+str(dataset_id)+'_plotdata.txt','w')
    output3.write('#RA, Dec, intFlx, intFlxErr, avgIntFlx, mjdTStart, Frequency, sepImgCentre, SideRealTime, Elevation, skymodel_flux, Image RMS, Image, runcatID \n')
    for x in range(len(source_data)):
        string='%s' % ','.join(str(val) for val in source_data[x])
        output3.write(string+'\n')
    output3.close()

source_data=tools.extract_data('ds'+str(dataset_id)+'_plotdata.txt')

unique_frequencies = np.unique([int(a[6]) for a in source_data])

print unique_frequencies

#for src in unique_sources:
#for x in sources:
#    jd = coords.julian_date(datetime.strptime(x[13].split('.')[0], '%Y-%m-%d %H:%M:%S').replace(tzinfo=pytz.utc), modified=True) # make time for the midpoint of the observation
#    if telescope == 'LOFAR' or telescope == 'MWA':
#        elevation,siderealtime = tools.sidetime(x[16], x[15], jd, telescope)
#       ### print elevation, siderealtime
#        
#        source_data.append([float(x[16]),float(x[15]),float(x[6]),float(x[7]),avgs[str(x[11])],jd,int(float(a[8])/1e6),coords.angsep(float(x[16]),float(x[15]),float(x[3]),float(x[2]))/3600.,siderealtime,elevation,src_assoc[str(x[11])], float(x[10]), x[14]]) # [RA, Dec, intFlx, intFlxErr, avgIntFlx, mjdTStart, Frequency, sepImgCentre, SideRealTime, Elevation, skymodel_flux, Image RMS, Image, runcatID]
#    else:
#        source_data.append([float(x[16]),float(x[15]),float(x[6]),float(x[7]),avgs[str(x[11])],coords.julian_date(datetime.strptime(x[13], '%Y-%m-%d %H:%M:%S.%f').replace(tzinfo=pytz.utc), modified=True),int(float(a[8])/1e6),coords.angsep(x[16],x[15],x[3],x[2])/3600.,-1,-1,src_assoc[str(x[11])], float(x[10]), x[14]]) # [RA, Dec, intFlx, intFlxErr, avgIntFlx, mjdTStart, Frequency, sepImgCentre, SideRealTime, Elevation, skymodel_flux, Image RMS, Image, runcatID]

plt.close()


# Setting up the plot
fig = plt.figure(1,figsize=(15,10))
nullfmt   = NullFormatter()         # no labels
fontP = FontProperties()
fontP.set_size('medium')

plt.subplots_adjust(wspace=0.0, hspace=0.0)

col = tools.make_colours(unique_frequencies)
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)

avg1 = np.average([float(x[2])/float(x[4]) for x in source_data])
stdev1 = np.std([float(x[2])/float(x[4]) for x in source_data])
avg2 = np.average([float(x[2])/float(x[10]) for x in source_data if float(x[10])!=-1])
stdev2 = np.std([float(x[2])/float(x[10]) for x in source_data if float(x[10])!=-1])


for i in range(len(unique_frequencies)):
    freq=unique_frequencies[i]
    
    ydata1 = [float(x[2])/float(x[4]) for x in source_data if int(x[6])==freq]
    ydata2 = [float(x[2])/float(x[10]) for x in source_data if int(x[6])==freq if float(x[10])!=-1]
    xdata1y1 = [float(x[2]) for x in source_data if int(x[6])==freq]
    xdata2y1 = [float(x[7]) for x in source_data if int(x[6])==freq]
    xdata1y2 = [float(x[2]) for x in source_data if int(x[6])==freq if float(x[10])!=-1]
    xdata2y2 = [float(x[7]) for x in source_data if int(x[6])==freq if float(x[10])!=-1]

    ax1.scatter(xdata1y1, ydata1,color=col[i], s=10.)
    ax2.scatter(xdata2y1, ydata1,color=col[i], s=10.)
    ax4.scatter(xdata1y2, ydata2,color=col[i], s=10.)
    ax5.scatter(xdata2y2, ydata2,color=col[i], s=10.)

ax1.legend(unique_frequencies, loc=2, prop=fontP)

ax3.hist([float(x[2])/float(x[4]) for x in source_data],bins=50,histtype='stepfilled', orientation='horizontal')
tmpData=[float(x[2])/float(x[10]) for x in source_data if float(x[10])!=-1]
if tmpData:
    ax6.hist(tmpData,bins=50,histtype='stepfilled', orientation='horizontal')

ax1.axhline(y=avg1, linewidth=2, color='k', linestyle='-')
ax2.axhline(y=avg1, linewidth=2, color='k', linestyle='-')
ax3.axhline(y=avg1, linewidth=2, color='k', linestyle='-')
ax4.axhline(y=avg2, linewidth=2, color='k', linestyle='-')
ax5.axhline(y=avg2, linewidth=2, color='k', linestyle='-')
ax6.axhline(y=avg2, linewidth=2, color='k', linestyle='-')

ax1.axhline(y=avg1+stdev1, linewidth=2, color='k', linestyle='--')
ax2.axhline(y=avg1+stdev1, linewidth=2, color='k', linestyle='--')
ax3.axhline(y=avg1+stdev1, linewidth=2, color='k', linestyle='--')
ax4.axhline(y=avg2+stdev2, linewidth=2, color='k', linestyle='--')
ax5.axhline(y=avg2+stdev2, linewidth=2, color='k', linestyle='--')
ax6.axhline(y=avg2+stdev2, linewidth=2, color='k', linestyle='--')

ax1.axhline(y=avg1-stdev1, linewidth=2, color='k', linestyle='--')
ax2.axhline(y=avg1-stdev1, linewidth=2, color='k', linestyle='--')
ax3.axhline(y=avg1-stdev1, linewidth=2, color='k', linestyle='--')
ax4.axhline(y=avg2-stdev2, linewidth=2, color='k', linestyle='--')
ax5.axhline(y=avg2-stdev2, linewidth=2, color='k', linestyle='--')
ax6.axhline(y=avg2-stdev2, linewidth=2, color='k', linestyle='--')

yvals = [float(x[2])/float(x[4]) for x in source_data]
xvals = [np.log10(float(x[2])) for x in source_data]
bins,yavg,ysig = tools.runningAvg(xvals,yvals)
ax1.plot(np.power(10.,bins),yavg,'r-', linewidth=2)

xvals = [float(x[7]) for x in source_data]
bins,yavg,ysig = tools.runningAvg(xvals,yvals)
ax2.plot(bins,yavg,'r-', linewidth=2)

yvals = [float(x[2])/float(x[10]) for x in source_data if float(x[10])!=-1]
xvals = [np.log10(float(x[2])) for x in source_data if float(x[10])!=-1]
if yvals:
    bins,yavg,ysig = tools.runningAvg(xvals,yvals)
    ax4.plot(np.power(10.,bins),yavg,'r-', linewidth=2)

    xvals = [float(x[7]) for x in source_data if float(x[10])!=-1]
    bins,yavg,ysig = tools.runningAvg(xvals,yvals)
    ax5.plot(bins,yavg,'r-', linewidth=2)




ax1.set_xscale('log')
ax4.set_xscale('log')
ax1.xaxis.set_major_formatter(nullfmt)
ax2.xaxis.set_major_formatter(nullfmt)
ax3.xaxis.set_major_formatter(nullfmt)
ax2.yaxis.set_major_formatter(nullfmt)
ax3.yaxis.set_major_formatter(nullfmt)
ax5.yaxis.set_major_formatter(nullfmt)
ax6.yaxis.set_major_formatter(nullfmt)
ax2.axes.yaxis.set_ticklabels([])
ax3.axes.yaxis.set_ticklabels([])
ax5.axes.yaxis.set_ticklabels([])
ax6.axes.yaxis.set_ticklabels([])
ax1.axes.xaxis.set_ticklabels([])
ax2.axes.xaxis.set_ticklabels([])
ax3.axes.xaxis.set_ticklabels([])
ax2.set_ylim( ax1.get_ylim() )
ax3.set_ylim( ax1.get_ylim() )
ax5.set_ylim( ax4.get_ylim() )
ax6.set_ylim( ax4.get_ylim() )
ax4.set_xlim( ax1.get_xlim() )
ax5.set_xlim( ax2.get_xlim() )
ax6.set_xlim( ax3.get_xlim() )

ax1.set_ylabel('Flux / Average Flux')
ax4.set_ylabel('Flux / Corrected Skymodel Flux')
ax4.set_xlabel('Flux (Jy)')
ax5.set_xlabel('Distance from centre of image (degrees)')
ax6.set_xlabel('Number')

plt.savefig('ds'+str(dataset_id)+'_flux_diagnostic1.png')
plt.close()

# Setting up the plot 2
fig = plt.figure(1,figsize=(15,10))
nullfmt   = NullFormatter()         # no labels
fontP = FontProperties()
fontP.set_size('medium')
col = tools.make_colours(unique_frequencies)

plt.subplots_adjust(wspace=0.0, hspace=0.0)

ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)

for i in range(len(unique_frequencies)):
    freq=unique_frequencies[i]
    
    ydata1 = [float(x[2])/float(x[4]) for x in source_data if int(x[6])==freq]
    ydata2 = [float(x[2])/float(x[10]) for x in source_data if int(x[6])==freq if float(x[10])!=-1]
    xdata1y1 = [float(x[11]) for x in source_data if int(x[6])==freq]
    xdata2y1 = [float(x[8]) for x in source_data if int(x[6])==freq]
    xdata3y1 = [float(x[9]) for x in source_data if int(x[6])==freq]
    xdata1y2 = [float(x[11]) for x in source_data if int(x[6])==freq if float(x[10])!=-1]
    xdata2y2 = [float(x[8]) for x in source_data if int(x[6])==freq if float(x[10])!=-1]
    xdata3y2 = [float(x[9]) for x in source_data if int(x[6])==freq if float(x[10])!=-1]

    ax1.scatter(xdata1y1, ydata1,color=col[i], s=10.)
    ax2.scatter(xdata2y1, ydata1,color=col[i], s=10.)
    ax3.scatter(xdata3y1, ydata1,color=col[i], s=10.)
    ax4.scatter(xdata1y2, ydata2,color=col[i], s=10.)
    ax5.scatter(xdata2y2, ydata2,color=col[i], s=10.)
    ax6.scatter(xdata3y2, ydata2,color=col[i], s=10.)


ax1.legend(unique_frequencies, loc=2, prop=fontP)


ax1.axhline(y=avg1, linewidth=2, color='k', linestyle='-')
ax2.axhline(y=avg1, linewidth=2, color='k', linestyle='-')
ax3.axhline(y=avg1, linewidth=2, color='k', linestyle='-')
ax4.axhline(y=avg2, linewidth=2, color='k', linestyle='-')
ax5.axhline(y=avg2, linewidth=2, color='k', linestyle='-')
ax6.axhline(y=avg2, linewidth=2, color='k', linestyle='-')

ax1.axhline(y=avg1+stdev1, linewidth=2, color='k', linestyle='--')
ax2.axhline(y=avg1+stdev1, linewidth=2, color='k', linestyle='--')
ax3.axhline(y=avg1+stdev1, linewidth=2, color='k', linestyle='--')
ax4.axhline(y=avg2+stdev2, linewidth=2, color='k', linestyle='--')
ax5.axhline(y=avg2+stdev2, linewidth=2, color='k', linestyle='--')
ax6.axhline(y=avg2+stdev2, linewidth=2, color='k', linestyle='--')

ax1.axhline(y=avg1-stdev1, linewidth=2, color='k', linestyle='--')
ax2.axhline(y=avg1-stdev1, linewidth=2, color='k', linestyle='--')
ax3.axhline(y=avg1-stdev1, linewidth=2, color='k', linestyle='--')
ax4.axhline(y=avg2-stdev2, linewidth=2, color='k', linestyle='--')
ax5.axhline(y=avg2-stdev2, linewidth=2, color='k', linestyle='--')
ax6.axhline(y=avg2-stdev2, linewidth=2, color='k', linestyle='--')

yvals = [float(x[2])/float(x[4]) for x in source_data]
xvals = [float(x[11]) for x in source_data]
bins,yavg,ysig = tools.runningAvg(xvals,yvals)
ax1.plot(bins,yavg,'r-', linewidth=2)

xvals = [float(x[8]) for x in source_data]
bins,yavg,ysig = tools.runningAvg(xvals,yvals)
ax2.plot(bins,yavg,'r-', linewidth=2)

xvals = [float(x[9]) for x in source_data]
bins,yavg,ysig = tools.runningAvg(xvals,yvals)
ax3.plot(bins,yavg,'r-', linewidth=2)

yvals = [float(x[2])/float(x[10]) for x in source_data if float(x[10])!=-1]
if yvals:
    xvals = [float(x[11]) for x in source_data if float(x[10])!=-1]
    bins,yavg,ysig = tools.runningAvg(xvals,yvals)
    ax4.plot(bins,yavg,'r-', linewidth=2)

    xvals = [float(x[8]) for x in source_data if float(x[10])!=-1]
    bins,yavg,ysig = tools.runningAvg(xvals,yvals)
    ax5.plot(bins,yavg,'r-', linewidth=2)

    xvals = [float(x[9]) for x in source_data if float(x[10])!=-1]
    bins,yavg,ysig = tools.runningAvg(xvals,yvals)
    ax6.plot(bins,yavg,'r-', linewidth=2)


    
ax1.xaxis.set_major_formatter(nullfmt)
ax2.xaxis.set_major_formatter(nullfmt)
ax3.xaxis.set_major_formatter(nullfmt)
ax2.yaxis.set_major_formatter(nullfmt)
ax3.yaxis.set_major_formatter(nullfmt)
ax5.yaxis.set_major_formatter(nullfmt)
ax6.yaxis.set_major_formatter(nullfmt)
ax2.axes.yaxis.set_ticklabels([])
ax3.axes.yaxis.set_ticklabels([])
ax5.axes.yaxis.set_ticklabels([])
ax6.axes.yaxis.set_ticklabels([])
ax1.axes.xaxis.set_ticklabels([])
ax2.axes.xaxis.set_ticklabels([])
ax3.axes.xaxis.set_ticklabels([])
ax2.set_ylim( ax1.get_ylim() )
ax3.set_ylim( ax1.get_ylim() )
ax5.set_ylim( ax4.get_ylim() )
ax6.set_ylim( ax4.get_ylim() )
ax4.set_xlim( ax1.get_xlim() )
ax5.set_xlim( ax2.get_xlim() )
ax6.set_xlim( ax3.get_xlim() )

ax1.set_ylabel('Flux / Average Flux')
ax4.set_ylabel('Flux / Corrected Skymodel Flux')
ax4.set_xlabel('Image RMS (Jy)')
ax5.set_xlabel('Sidereal time (hours)')
ax6.set_xlabel('Elevation (degrees)')

plt.savefig('ds'+str(dataset_id)+'_flux_diagnostic2.png')
plt.close()









