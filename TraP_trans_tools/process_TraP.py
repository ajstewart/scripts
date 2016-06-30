import format_TraP_data
import plotting_tools
import generic_tools
import numpy as np
import sys
import os

# Obtain input parameters from the command line
if len(sys.argv) != 11:
    print 'python process_TraP.py <database> <username> <password> <dataset_id> <release> <host> <port> <sigma1> <sigma2> <lightcurves>'
    exit()
database = sys.argv[1]
username = sys.argv[2]
password = sys.argv[3]
dataset_id = str(sys.argv[4])
release = str(sys.argv[5])
host = str(sys.argv[6])
port = int(sys.argv[7])
sigma1 = float(sys.argv[8])
sigma2 = float(sys.argv[9])
lightcurves = sys.argv[10]

# get TraP data from the database and sort it into the required array which is then loaded
if not os.path.isfile('ds'+str(dataset_id)+'_trans_data.txt'):
    format_TraP_data.format_data(database,dataset_id,release,host,port,username,password,lightcurves)
trans_data, transdata2=generic_tools.extract_data('ds'+str(dataset_id)+'_trans_data.txt')
# make first array for the scatter_hist plot: [log10(eta_nu), log10(V_nu), nu]
data=[[trans_data[n][0],np.log10(float(trans_data[n][1])),np.log10(float(trans_data[n][2])),trans_data[n][6], trans_data[n][-2]] for n in range(len(trans_data)) if float(trans_data[n][1]) > 0 if float(trans_data[n][2]) > 0 if trans_data[n][-5]=='2']

# print out the transients that TraP automatically found
print 'Identified Transient Candidates (no margin)'
transcandmask=(transdata2[:,-5]!='2') & (transdata2[:,-3].astype(float)>=transdata2[:,-2].astype(float)) & (transdata2[:,-4].astype(float)<transdata2[:,-2].astype(float))
transcands=transdata2[transcandmask,:]
print np.sort(transcands[:,0])
print "Total number of Transient Candidates (no margin) = {0}".format(len(transcands[:,0]))
# print np.sort(list(set([int(x[0]) for x in trans_data if x[-4]!='2' if float(x[-2])>=float(x[-1]) if float(x[-3])<float(x[-1])])))
np.savetxt('transient_candidates.txt', transcands, fmt='%s', delimiter=',', header='Runcat_id, eta_nu, V_nu, max_flux, mean_flux, fluxrat, freq, dpts, RA, Dec, date, trans_type, max_rms_sigma, min_rms_sigma, detection_threshold, sig_to_noise')

print 'Identified Transients (no margin)'
transmask=(transdata2[:,-5]!='2') & (transdata2[:,-4].astype(float)>=transdata2[:,-2].astype(float))
trans=transdata2[transmask,:]
print np.sort(trans[:,0])
print "Total number of Transients (no margin) = {0}".format(len(trans[:,0]))
# print np.sort(list(set([int(x[0]) for x in trans_data if x[-4]!='2' if float(x[-3])>=float(x[-1])])))
np.savetxt('transients.txt', trans, fmt='%s', delimiter=',', header='Runcat_id, eta_nu, V_nu, max_flux, mean_flux, fluxrat, freq, dpts, RA, Dec, date, trans_type, max_rms_sigma, min_rms_sigma, detection_threshold, sig_to_noise')

# Find the thresholds for a given sigma (in log space)
sigcutx,paramx,range_x = generic_tools.get_sigcut([a[1] for a in data],sigma1)
sigcuty,paramy,range_y = generic_tools.get_sigcut([a[2] for a in data],sigma2)
if sigma1 == 0:
    sigcutx=0
if sigma2 == 0:
    sigcuty=0
print(r'Gaussian Fit $\eta$: '+str(round(10.**paramx[0],2))+'(+'+str(round((10.**(paramx[0]+paramx[1])-10.**paramx[0]),2))+' '+str(round((10.**(paramx[0]-paramx[1])-10.**paramx[0]),2))+')')
print(r'Gaussian Fit $V$: '+str(round(10.**paramy[0],2))+'(+'+str(round((10.**(paramy[0]+paramy[1])-10.**paramy[0]),2))+' '+str(round((10.**(paramy[0]-paramy[1])-10.**paramy[0]),2))+')')
print 'Eta_nu threshold='+str(10.**sigcutx)+', V_nu threshold='+str(10.**sigcuty)

# Get the different frequencies in the dataset
frequencies = generic_tools.get_frequencies(data)

# Create the scatter_hist plot
IdTrans, IdTrans2, IdTrans3 = plotting_tools.create_scatter_hist(data,sigcutx,sigcuty,paramx,paramy,range_x,range_y,dataset_id,frequencies)

print 'Identified variables:'
if len(IdTrans)>0:
    print np.sort(IdTrans[:,0])
    variablemask=np.in1d(transdata2[:,0], IdTrans[:,0])
    variablearray=transdata2[variablemask,:]
else:
    variablearray=[]
    print variablearray
print "Total number of variables = {0}".format(len(variablearray))
np.savetxt('variables.txt', variablearray, fmt='%s', delimiter=',', header='Runcat_id, eta_nu, V_nu, max_flux, mean_flux, fluxrat, freq, dpts, RA, Dec, date, trans_type, max_rms_sigma, min_rms_sigma, detection_threshold, sig_to_noise')
# np.savetxt('variables.txt', variablearray[variablearray[:,0].argsort()], fmt='%s', delimiter=',', header='Runcat_id, eta_nu, V_nu, flux, fluxrat, freq, dpts, RA, Dec, trans_type, max_rms_sigma, min_rms_sigma, detection_threshold')
# f=open('variables.txt', 'w')
# for i in IdTrans:
#     f.write(i+'\n')
# f.close()

print 'Identified variable candidates (eta):'
if len(IdTrans2)>0:
    print np.sort(IdTrans2[:,0])
    variableetamask=np.in1d(transdata2[:,0], IdTrans2[:,0])
    variableetaarray=transdata2[variableetamask,:]
else:
    variableetaarray=[]
    print variableetaarray
print "Total number of variable candidates (eta) = {0}".format(len(variableetaarray[:,0]))
np.savetxt('variables_eta.txt', variableetaarray[variableetaarray[:,0].argsort()], fmt='%s', delimiter=',', header='Runcat_id, eta_nu, V_nu, max_flux, mean_flux, fluxrat, freq, dpts, RA, Dec, date, trans_type, max_rms_sigma, min_rms_sigma, detection_threshold, sig_to_noise')
# print np.sort(list(set(IdTrans2)))
# f=open('variables_eta.txt', 'w')
# for i in IdTrans2:
#     f.write(i+'\n')
# f.close()
#
print 'Identified variable candidates (V):'
if len(IdTrans3)>0:
    print np.sort(IdTrans3[:,0])
    variableVmask=np.in1d(transdata2[:,0], IdTrans3[:,0])
    variableVarray=transdata2[variableVmask,:]
else:
    variableVarray=[]
    print variableVarray
print "Total number of variable candidates (V) = {0}".format(len(variableVarray[:,0]))
np.savetxt('variables_V.txt', variableVarray[variableVarray[:,0].argsort()], fmt='%s', delimiter=',', header='Runcat_id, eta_nu, V_nu, max_flux, mean_flux, fluxrat, freq, dpts, date, RA, Dec, trans_type, max_rms_sigma, min_rms_sigma, detection_threshold')
# print np.sort(list(set(IdTrans3)))
# f=open('variables_V.txt', 'w')
# for i in IdTrans3:
#     f.write(i+'\n')
# f.close()

# make second array for the diagnostic plot: [eta_nu, V_nu, maxflx_nu, flxrat_nu, nu, trans_type]
data2=[[trans_data[n][0],float(trans_data[n][1]),float(trans_data[n][2]),float(trans_data[n][3]),float(trans_data[n][5]),trans_data[n][6], trans_data[n][-2]] for n in range(len(trans_data)) if float(trans_data[n][1]) > 0 if float(trans_data[n][2]) > 0 if trans_data[n][-5]=='2'] 

# Create the diagnostic plot
plotting_tools.create_diagnostic(data2,sigcutx,sigcuty,frequencies,dataset_id)
