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

work_dir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/scratch'
fname = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
db = gct.GCT()
db.read(fname)

##load query result - gctx file
rslt = gct.GCT()
#if specific result directory is specified, use that - otherwise get gctx from working dir
rslt.read('/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/jan04/dose_analysis_tool.1357311906761/result_WTESLM.COMBINED_n85x398050.gctx')

### signature consistancy graph ### 
qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qDose = db.get_column_meta('pert_dose')
probeIDs = db.get_row_meta('id')
#set null distirbution of z-scores (currently normal)
ESmat = db.matrix
#calculate null distribution
mu, sigma = 0, 1
s = numpy.random.normal(mu, sigma, len(ESmat[:,1]))
s.sort()

pertSet = set(qPert)
#pertSet = set(qPert[10:20])
ProbePercMtrx = numpy.zeros(shape=(978,len(pertSet)))
for i1, pert in enumerate(pertSet):
	iP = _all_indices(pert, qPert) #index of doses on plate
	if len(iP) < 2:
		print pert + ' has only one instance'
		# continue
	uDose = [qDose[i] for i in iP]
	fDose = [float(x) for x in uDose] #convert strings to float
	aDose = numpy.asarray(fDose) #convert to numpy array
	iD = aDose.argsort() #local ordering
	sDose = [fDose[j] for j in iD] #sort local doses
	iPo =  [iP[i] for i in iD] #ordered index
	#sMat = ESmat[:,iPo]
	#sMat.sort(axis=0)
	#mongo query for each unique pertID
	qStr = qPertID[iPo[0]] #set pertID
	if len(qStr) >= 13:
		qStr = qStr[0:13] #shorten qPertID
	CM = mutil.CMapMongo()
	#cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True})
	# edge50Lst = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True,'up50_lm':True,'dn50_lm':True,'cell_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
	edge50Lst = CM.find({'pert_id':{'$regex':qStr},'cell_id':'PC3'},{'sig_id':True,'up50_lm':True,'dn50_lm':True,'cell_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
	nInstances = len(edge50Lst) #number of instances in db
	#count number of times a probe is in the top/bottom 50 genes of an instance
	upProbeCnts = [0] * len(probeIDs)
	dnProbeCnts = [0] * len(probeIDs)
	for j,inst in enumerate(edge50Lst):
		up50 = edge50Lst[j]['up50_lm']
		dn50 = edge50Lst[j]['dn50_lm']
		#loop through every gene in the top and bottom list - where does it live on the rank list?
		for prb in up50:
			if prb in probeIDs:
				iPrb = probeIDs.index(prb)
				upProbeCnts[iPrb] = upProbeCnts[iPrb] +1
		for prb in dn50:
			if prb in probeIDs:
				iPrb = probeIDs.index(prb)
				dnProbeCnts[iPrb] = dnProbeCnts[iPrb] +1
	upCnts = numpy.asarray(upProbeCnts)
	upPerct = upCnts/float(nInstances)
	upPerct.sort()
	ProbePercMtrx[:,i1] = upPerct

#loop loop though each unique compound and plot
fig = plt.figure(figsize=(8.0, 2.0))
ax = fig.add_subplot(111)
for i, v in enumerate(pertSet):
	ax.plot(range(len(upPerct)),ProbePercMtrx[:,i],color='r')
ax.plot(range(len(upPerct)),ProbePercMtrx[:,21],'o-',color='b')
ax.set_xlabel('ordered probes')
ax.set_ylabel('percent of instance edges')
ax.set_title('frequency of probes occuring in edge list for a given compound')
plt.show()



#convert to percentage of total instances
fig = plt.figure(figsize=(8.0, 2.0))
ax = fig.add_subplot(111)
# # the histogram of the data
# n, bins, patches = ax.hist(upPerct, 50, facecolor='green', alpha=0.75)
# ax.bar(range(len(upPerct)),upPerct)
ax.plot(range(len(upPerct)),upPerct)
plt.show()
