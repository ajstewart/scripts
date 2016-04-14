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
import ConfigParser
import optparse
import getpass
import socket

###################### INITIAL SETUP STEPS ######################

user=getpass.getuser()
hostname=socket.gethostname()

usage = "usage: python %prog [options] "
description="A script that uses TraP outputs analyse the quality of images and outputs the recommended QC settings for TraP. \
Script updated for TraP Release 2 databases. If using TraP Release 1.1 please use older version."
vers="2.0"

parser = optparse.OptionParser(usage=usage, version="%prog v{0}".format(vers), description=description)
parser.add_option("--config", action="store", type="string", dest="config", default="pipeline.cfg", help="Optionally define the pipeline.cfg file to use from your TraP set up to\
quickly pass database details [default: %default]")
parser.add_option("-d", "--database", action="store", type="string", dest="database", default=user,help="The name of the TraP database you are using [default: %default]")
parser.add_option("-u", "--username", action="store", type="string", dest="username", default=user, help="Your username for the database [default: %default]")
parser.add_option("-p", "--password", action="store", type="string", dest="password", default=user, help="The password for the database [default: %default]")
parser.add_option("-H", "--host", action="store", type="string", dest="host", default=hostname, help="The name of the machine hosting the databases [default: %default]")
parser.add_option("-o", "--port", action="store", type="string", dest="port", default="5432", help="The port number for the machine, typically 5432 for postgresql and 52000 for monetdb [default: %default]")
parser.add_option("-t", "--databasetype", action="store", type="choice", choices=["postgresql", "monetdb"], dest="databasetype", default="postgresql", help="postgresql or monetdb [default: %default]")
parser.add_option("-i", "--datasetid", action="store", type="string", dest="datasetid", default="1", help="The dataset containing all the images rejected [default: %default]")
parser.add_option("-s", "--sigma", action="store", type="float", dest="sigma", default=2.0, help="The sigma clipping to be used for the RMS highbound [default: %default]")
parser.add_option("-f", "--plotfreqs", action="store_true", dest="plotfreqs", default=False, help="Use option to plot all QC plots for individual frequencies [default: %default]")
parser.add_option("-x", "--datasetid2", action="store", type="string", dest="datasetid2", default="N", help="The dataset containing all the extracted sources for the images\
 - if you do not have this yet, leave as the default value N [default: %default]")

(options, args) = parser.parse_args()

configfile=options.config

if os.path.isfile(configfile):
    config = ConfigParser.ConfigParser()
    config.read(configfile)
    database = config.get("database", "database")
    username = config.get("database", "user")
    password = config.get("database", "password")
    host = config.get("database", "host")
    port = config.get("database", "port")
    databaseType = config.get("database", "engine")
    ports={"postgresql":"5432", "monetdb":"52000"}
    if port =="":
        port=ports[databaseType]
else:
    database = options.database
    username = options.username
    password = options.password
    host = options.host
    port = options.port
    databaseType = options.databasetype

# if len(sys.argv) != 11:
#     print 'python TraP_QC_diagnostics.py <database> <username> <password> <host> <port> <databaseType> <dataset_id> <sigma> <plt_freqs> <database_id2>'
#     exit()
dataset_id = options.datasetid
sigma = options.sigma
plt_freqs = options.plotfreqs
dataset_id2 = options.datasetid2
# A-Team positions
CasA=[350.866417,58.811778]
CygA=[299.868153,40.733916]
VirA=[187.705930,12.391123]

min_sep=0. # The absolute minimum allowed separation from the A-Team source, set to zero to enable the code to work independently

###################### MAIN SCRIPT ######################

# Extracting data from the TraP database into a text file
if not os.path.exists('ds_'+dataset_id+'_images.csv'):
    try:
        tools.get_data(database, username, password, host, port, databaseType, dataset_id, dataset_id2)
    except:
        print "Cannot connect to database with settings:\n"
        print "database = {0}".format(database)
        print "username = {0}".format(username)
        print "password = {0}".format(password)
        print "host = {0}".format(host)
        print "port = {0}".format(port)
        print "database_type = {0}".format(databaseType)
        print "dataset_id = {0}".format(dataset_id)
        print "\nPlease check database settings and/or dataset id."
        sys.exit()

# Extract relevant data from dataset text file
image_info, frequencies, plt_ratios = tools.extract_data(dataset_id, CasA, CygA, VirA)

freq='all'
# RMS noise properties
noise_avg_log, noise_scatter_log, noise_threshold_log = tools.fit_hist([np.log10(image_info[n][4]*1e3) for n in range(len(image_info))], sigma, r'Observed RMS (mJy)', 'ds'+dataset_id+'_rms', freq)
noise_avg=10.**(noise_avg_log)
noise_max=10.**(noise_avg_log+noise_scatter_log)-10.**(noise_avg_log)
noise_min=10.**(noise_avg_log)-10.**(noise_avg_log-noise_scatter_log)
print 'Average RMS Noise in images (1 sigma range, frequency='+str(freq)+' MHz): '+str(round(noise_avg,1))+' (+'+str(round(noise_max,1))+',-'+str(round(noise_min,1))+') mJy'
if plt_ratios:
    # RMS/Theoretical limit for TraP
    ratio_avg_log, ratio_scatter_log, ratio_threshold_log = tools.fit_hist([np.log10(image_info[n][6]) for n in range(len(image_info))], sigma, r'Observed RMS / Theoretical Noise', 'ds'+dataset_id+'_ratio', freq)
    ratio_avg=10.**(ratio_avg_log)
    ratio_threshold = round((10.**ratio_threshold_log),1)
    ratio_threshold2 = round((10.**((ratio_avg_log - ratio_threshold_log)+ratio_avg_log)),1)
    print 'Average RMS/Theoretical in images (frequency='+str(freq)+' MHz): '+str(round(ratio_avg,1))
    print '######## Recommended TraP high_bound threshold: '+str(ratio_threshold)
    print '######## Recommended TraP low_bound threshold: '+str(ratio_threshold2)
    ratio_avg_log, ratio_scatter_log, ratio_threshold_log = tools.fit_hist([np.log10(image_info[n][-1]) for n in range(len(image_info))], sigma, r'Observed RMS / Confusion Noise', 'ds'+dataset_id+'_confratio', freq)

else:
    ratio_threshold=10000000.
    ratio_threshold2=0.

tools.plotfig_scatter(image_info, 7, 4, 'Ellipticity (Bmaj/Bmin)', 'RMS (Jy)', 'ds'+dataset_id+'_theoretical_ellipticity_'+str(freq)+'MHz')

if plt_freqs:
    for freq in frequencies:
        # RMS noise properties
        noise_avg_log_tmp, noise_scatter_log_tmp, noise_threshold_log_tmp = tools.fit_hist([np.log10(image_info[n][4]) for n in range(len(image_info)) if image_info[n][3]==freq], sigma, r'Observed RMS (Jy)', 'ds'+dataset_id+'_rms', freq)
        noise_avg_tmp=10.**(noise_avg_log_tmp)
        noise_max_tmp=10.**(noise_avg_log_tmp+noise_scatter_log_tmp)-10.**(noise_avg_log_tmp)
        noise_min_tmp=10.**(noise_avg_log_tmp)-10.**(noise_avg_log_tmp-noise_scatter_log_tmp)
        print 'Average RMS Noise in images (1 sigma range, frequency='+str(freq)+' MHz): '+str(round(noise_avg_tmp*1e3,1))+' (+'+str(round(noise_max_tmp*1e3,1))+',-'+str(round(noise_min_tmp*1e3,1))+') mJy'
        if plt_ratios:
            # RMS/Theoretical limit for TraP
            ratio_avg_log_tmp, ratio_scatter_log_tmp, ratio_threshold_log_tmp = tools.fit_hist([np.log10(image_info[n][6]) for n in range(len(image_info)) if image_info[n][3]==freq], sigma, r'Observed RMS / Theoretical Noise', 'ds'+dataset_id+'_ratio', freq)
            ratio_avg_tmp=10.**(ratio_avg_log_tmp)
            print 'Average RMS/Theoretical in images (frequency='+str(freq)+' MHz): '+str(round(ratio_avg_tmp,1))

# Calculate restoring beam threshold using a simple clipping using the average and rms scatter
rms2=[image_info[x][7] for x in range(len(image_info))]
avg_rms2=(sum(rms2)/len(rms2))
rms_rms2=math.sqrt((sum(n*n-(avg_rms2*avg_rms2) for n in rms2))/len(rms2))
ellipticity_threshold=round(avg_rms2+sigma*rms_rms2,2)
print 'Average ellipticity: '+str(round(avg_rms2,2))+' +/- '+str(round(rms_rms2,2))
print '######## Recommended TraP ellipticity threshold: '+str(ellipticity_threshold)

image_info_clip1 = [image_info[n] for n in range(len(image_info)) if (image_info[n][6] < ratio_threshold) and (image_info[n][7] < ellipticity_threshold)]
images_rejected_clip1 = ["#Images rejected on basis of {0} sigma RMS ratio cut ({1}) and recommended ellipticity cut ({2})".format(sigma, ratio_threshold, ellipticity_threshold)]+sorted([image_info[n][-4] for n in range(len(image_info)) if image_info[n][6] > ratio_threshold or image_info[n][7] > ellipticity_threshold])

for ateam in ['CasA', 'CygA', 'VirA']:
    if ateam == 'CasA':
        a=8
    elif ateam == 'CygA':
        a=9
    elif ateam == 'VirA':
        a=10
    min_sep = tools.plotfig_ATeam(image_info_clip1, a, 4, 'Separation (degrees)', 'RMS (Jy/beam)', 'ds'+dataset_id+'_'+ateam+'_clipped',min_sep,noise_avg)
print '######## Recommended TraP min seperation from A-Team sources (degrees): '+str(min_sep)

image_info_clip2 = [image_info_clip1[n] for n in range(len(image_info_clip1)) if image_info_clip1[n][8] > min_sep if image_info_clip1[n][9] > min_sep if image_info_clip1[n][10] > min_sep]
images_rejected_clip2 = ["#Images rejected on basis of recommended distance from A-team sources ({0} deg)".format(min_sep),]+sorted([image_info_clip1[n][-4] for n in range(len(image_info_clip1)) if image_info_clip1[n][8] < min_sep if image_info_clip1[n][9] < min_sep if image_info_clip1[n][10] < min_sep])
tools.plotfig_scatter(image_info_clip2, 7, 4, 'Ellipticity (Bmaj/Bmin)', 'RMS (Jy)', 'ds'+dataset_id+'_theoretical_ellipticity_allMHz_Final')

print 'Total images: '+str(len(image_info))+' After first clip: '+str(len(image_info_clip1))+' ('+str(int(round(100.*(float(len(image_info_clip1))/float(len(image_info))),0)))+'%) After second clip: '+str(len(image_info_clip2))+' ('+str(int(round(100.*(float(len(image_info_clip2))/float(len(image_info))),0)))+'%)'

np.savetxt("ds"+str(dataset_id)+"_image_info.csv", image_info, fmt='%s', delimiter=",")
np.savetxt("ds"+str(dataset_id)+"_images_rejected_clip1.txt", images_rejected_clip1, fmt='%s')
np.savetxt("ds"+str(dataset_id)+"_images_rejected_clip2.txt", images_rejected_clip2, fmt='%s')

###################### IMAGES AVAILABLE? ######################

if dataset_id2=='N':
    print 'No sources available'
    exit()

avg_flxrat=[]
sources = tools.extr_src_data(dataset_id2)
sky_data={}
flux_data=[]
for img in [x for x in image_info_clip2]:
    srcs=[x for x in sources if x[0] == img[11]]
    skymodel=str(int(round(float(img[1]),0)))+'_'+str(int(round(float(img[0]),0)))+'.sky'
    skymodel_key=str(int(round(float(img[1]),0)))+str(int(round(float(img[0]))))
    search_radius=5.
    flux_limit=0.
    if not os.path.isfile(skymodel):
        os.system('gsm.py '+skymodel+' '+str(img[1])+' '+str(img[0])+' '+str(search_radius)+' '+str(flux_limit))
    if skymodel_key not in sky_data.keys():
        vlss_sources=[]
        vlss_data=open(skymodel, 'r')
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
    vlss=tools.extract_sky(sky_data[skymodel_key],img[3],frq)
    tmp= tools.source_assoc(vlss,srcs,img[-1])
    for t in range(len(tmp)):
        flux_data.append([img[2],img[3],tmp[t][0], tmp[t][1]])
    flxrat, flxrms = tools.find_avg_int_flx_rat(vlss,srcs,img[-1])
    avg_flxrat.append([img[4],flxrat,img[3]])
freq='all'
flx_avg_log_tmp, flx_scatter_log_tmp, flx_threshold_log_tmp = tools.fit_hist([x[1] for x in avg_flxrat if x[1]!=0.], 0.0, r'Average (Flux / Corrected Skymodel Flux)', 'ds'+dataset_id+'_flux', freq)
print 'Average Flux Ratio: '+str(flx_avg_log_tmp)+' +/-'+str(flx_scatter_log_tmp)
tools.plotfig_scatter([x for x in avg_flxrat if x[1] != 0.], 0, 1, 'RMS (Jy)', 'Average(Flux / Corrected Skymodel Flux)', 'ds'+dataset_id+'_flux_'+str(freq)+'MHz_Final')
if plt_freqs:
    for freq in frequencies:
        tools.plotfig_scatter([x for x in avg_flxrat if x[2] == freq if x[1]!=0.], 0, 1, 'RMS (Jy)', 'Average(Flux / Corrected Skymodel Flux)', 'ds'+dataset_id+'_flux_'+str(freq)+'MHz_Final')
        flx_avg_log_tmp, flx_scatter_log_tmp, flx_threshold_log_tmp = tools.fit_hist([x[1] for x in avg_flxrat if x[1]!=0. if x[2]==freq], 0.0, r'Average (Flux / Corrected Skymodel Flux)', 'ds'+dataset_id+'_flux', freq)
        print 'Average Flux Ratio ('+str(freq)+' MHz): '+str(flx_avg_log_tmp)+' +/-'+str(flx_scatter_log_tmp)

np.savetxt("ds"+str(dataset_id)+"_flux_info.csv", flux_data, fmt='%s', delimiter=",")
