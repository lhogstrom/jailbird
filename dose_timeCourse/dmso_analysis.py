#! /usr/bin/env python
'''
to evaluate the null distribution of z-scores in DMSO in specific cell lines
'''

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


# qStr ='BRD-A19037878-001-04-9' #trichostatin-a
qStr = 'DMSO'
cell1 = 'PC3'
CM = mutil.CMapMongo()
dmsoLst = CM.find({'pert_id':{'$regex':qStr},'cell_id':cell1},{'sig_id':True})
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
probeIDs = affogato.get_row_meta('id')

### QQ plot ###
# qStr = 'DMSO'
edge50Lst = CM.find({'pert_id':{'$regex':qStr},'cell_id':cell1},{'sig_id':True,'up50_lm':True,'dn50_lm':True,'cell_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
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
for d in range(10):
#count probe enrichment and plot
		cmpd1 = qStr
		# dose1 = qDose[d]
		zLst = affogato.matrix[:,d]
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
		#plt.legdend([b1], ['probe in 100% of ' + str(nInstances) + 'instances' ], numpoints=1)
		# ax.set_title(pert + ' dose = ' + dose1)
		fname = pert + '_' + dose1 + 'um_connection_qq.png'
		outf = os.path.join(work_dir,fname)
		print 'saving graph for ' + pert 
		# plt.savefig(outf, bbox_inches=0)


### are there specific genes that get up/down regulated in DMSO (dimethyl sulfoxide)
upPercMkrs = numpy.divide(upProbeCnts,float(nInstances))
plt.hist(upPercMkrs,50)
plt.xlabel('relative occurance of a probe in the up50_lm list')
plt.ylabel('freq - number of probes')
plt.show()
