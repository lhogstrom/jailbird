import os
import jinja2
import argparse
import numpy
import scipy.stats as stats
import random
import matplotlib.pyplot as plt
import cmap
import cmap.analytics.sc as sc
import cmap.io.gct as gct
import cmap.util.tool_ops as tool_ops
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import cmap.analytics.signature_strength as ss
import glob
import subprocess
import time
import re
import operator
import multiprocessing

def _all_indices(value, qlist):
	'''
	return all indices where the input string matches the item in the list
	'''
	indices = []
	indx = -1
	while True:
		try:
			indx = qlist.index(value, indx+1)
			indices.append(indx)
		except ValueError:
			break
	# indices = [i for i, x in enumerate(qlist) if x == value]
	return indices

def _log_gen(x):
	'''
	generator function to build log series.  x controls the size of the log
	step in the series. x=1 yields full log steps, x=0.5 yields
	half log steps, etc.
	'''
	a = 1
	while True:
		yield a
		a *= 10**x

work_dir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/scratch'
fname = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
db = gct.GCT()
db.read(fname)

qStr = 'DMSO'
CM = mutil.CMapMongo()
dmsoLst = CM.find({'pert_id':{'$regex':qStr},'cell_id':'PC3'},{'sig_id':True})
nDmsos = len(dmsoLst)
# #save dmso list
# dmsoOut = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/scratch/dmsoPC3_list.grp'
# with open(dmsoOut,'w') as f:
# 	for d in dmsoLst:
# 		f.write(d + '\n')
# load in affogato matrix data
affogato = gct.GCT()
fAff = '/xchip/cogs/data/build/affogato/affogato_r1_score_n398050x22268.gctx'
LMrid = db.get_rids()
affogato.read(fAff,cid=dmsoLst,rid=LMrid)
mtrx = affogato.matrix


# loop through and test the null distribution of expression-template correlations with different dose numbers
i=0
m = 20 #number of theoritical doses to test
test_doses=[4,7,10,13]
m=len(test_doses)
LinCorrMtrx = numpy.zeros(shape=(978,m))
for n in test_doses:
	#pull randomly from the DMSO
	# n = 4 #number of cids to pull
	iN = numpy.random.randint(0,nDmsos,n)
	nMtrx = mtrx[:,iN]

	# build prototype curves if there is more than one dose
	linear = numpy.linspace(1,10,n)
	log_gen = _log_gen(1)
	log_curve = [log_gen.next() for x in range(n)]
	log_gen = _log_gen(.5)
	half_log_curve = [log_gen.next() for x in range(n)]
	log_gen = _log_gen(.25)
	quarter_log_curve = [log_gen.next() for x in range(n)]
	curves = numpy.array([linear,log_curve,
						  half_log_curve,quarter_log_curve])
	num_probes = nMtrx.shape[0]
	cc = numpy.corrcoef(nMtrx,curves)

	# grab the correlation values for all the probes against prototype curves
	linear_probe_corrs = cc[0:num_probes,num_probes]
	log_probe_corrs = cc[0:num_probes,num_probes + 1]
	half_log_probe_corrs = cc[0:num_probes,num_probes + 2]
	quarter_log_probe_corrs = cc[0:num_probes,num_probes + 3]

	#save correlation data for all compounds to a matrix
	LinCorrMtrx[:,i] = linear_probe_corrs 
	i = i+1

plt.interactive(True)
for j in range(LinCorrMtrx.shape[1]):
	col1 = LinCorrMtrx[:,j]
	col1.sort()
	plt.plot(col1,'.')
plt.legend(['4 doses', '7 doses', '10 doses', '13 doses'],numpoints=1, loc=4)
plt.xlabel('probes sorted by correlation with template')
plt.ylabel('corr with linear dose template')
plt.show()






### dose analysis 
# qq plot
cmpd1 = qPert[d]
dose1 = qDose[d]
zLst = db.matrix[:,d]
iLst = zLst.argsort() #sort z-scores and save index
sLst = zLst[iLst]
sUpProbeCnts = [upProbeCnts[l] for l in iLst] #sort probe counts acording to z-score
sDnProbeCnts = [dnProbeCnts[l] for l in iLst]
#mkrs = numpy.sqrt(sprobeCnts) # non linear scaling of marker points
sUpProbeCnts = [float(l) for l in sUpProbeCnts] #convert to float
sDnProbeCnts = [float(l) for l in sDnProbeCnts] #convert to float
# upPercMkrs = numpy.divide(sUpProbeCnts,max(sUpProbeCnts)) #divide by max count to make for relative frequency
# dnPercMkrs = numpy.divide(sDnProbeCnts,max(sDnProbeCnts))
upPercMkrs = numpy.divide(sUpProbeCnts,nInstances) #divide by total instances to make for relative frequency
dnPercMkrs = numpy.divide(sDnProbeCnts,nInstances)
upMkrs = numpy.multiply(upPercMkrs,100)
dnMkrs = numpy.multiply(dnPercMkrs,100)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(s,s,'b')
for j,sl in enumerate(sLst):
	ax.plot(s[j],sl,'r.',markersize=upMkrs[j],alpha=.25)
	ax.plot(s[j],sl,'b.',markersize=dnMkrs[j],alpha=.25)
ax.set_ylabel('observed z-score')
ax.set_xlabel('expected z-score')
# #set legend based on the number of
r1 = ax.plot(0,0,'r.',markersize=100,alpha=.25)
b1 = ax.plot(0,0,'b.',markersize=100,alpha=.25)
legStrUp = 'probe in 100% of ' + str(nInstances) + ' UP instances'
legStrDn = 'probe in 100% of ' + str(nInstances) + ' DN instances'
plt.legend([r1, b1], [legStrUp, legStrDn], numpoints=1, loc=4)





fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
r1 = ax2.plot(0,0,'r.',markersize=100,alpha=.25)
b1 = ax2.plot(0,0,'b.',markersize=100,alpha=.25)
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(s[j],sl,'r.',markersize=upMkrs[j],alpha=.25)
ax.plot(s[j],sl,'b.',markersize=dnMkrs[j],alpha=.25)
plt.legend([r1, b1], [legStrUp, legStrDn], numpoints=1, loc=4)




ax.set_ylabel('observed z-score')
ax.set_xlabel('expected z-score')
# #set legend based on the number of
r1 = ax.plot(0,0,'r.',markersize=100,alpha=.25)
b1 = ax.plot(0,0,'b.',markersize=100,alpha=.25)
legStrUp = 'probe in 100% of ' + str(nInstances) + ' UP instances'
legStrDn = 'probe in 100% of ' + str(nInstances) + ' DN instances'
plt.legend([r1, b1], [legStrUp, legStrDn], numpoints=1, loc=4)



#load in unique pert Ids with pert Desc
pertIDs = db.get_column_meta('pert_id')
pertDescs = db.get_column_meta('pert_desc')
upertIDs = set(pertIDs)
upertDescs = []
for pID in upertIDs:
	i = pertIDs.index(pID)
	upertDescs.append(pertDescs[i])
	
for x in upertIDs:
	print x[:13]
for x in upertDescs:
	print x 

