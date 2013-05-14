#!/bin/py
'''
analyze dose data with prism
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
import re
import scipy.stats
import pylab

#return all indices where the input string matches the item in the list
def _all_indices(value, qlist):
	'''
	input: 1) string and 2) list - return all indices where the input string matches the item in the list
	'''
	indices = []
	indx = -1
	while True:
		try:
			indx = qlist.index(value, indx+1)
			indices.append(indx)
		except ValueError:
			break
	return indices

#load in prism plates
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism'

#run with all cell lines/ time points
cellLst = ['PC3', 'A375', 'MCF7']
timeLst = ['6H', '24H']
#run with just one cell line/ time point
# cellLst = ['PC3']
# timeLst = ['6H']

for cell in cellLst:
	for tim in timeLst:
		### make SC plots
		cellLine = cell
		timeP = tim
		refControl = 'pc' #use pc vs vc controled data
		gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/by_pert_id_pert_dose/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		gctfile = gctfile[0]
		work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/%s_%s_%s' % (cell,timeP,refControl)
		if not os.path.exists(work_dir):
			os.mkdir(work_dir)
		db = gct.GCT() #make a gct object
		db.read(gctfile)
		### copy stuff from query tool 
		# instantiate a progress object
		prog = progress.DeterminateProgressBar('Dose Analysis')
		# make an SC object from the given gctx file
		sco = sc.SC()
		sco.add_sc_from_gctx_meta(gctfile, verbose=False)
		sco.set_thresh_by_specificity(0.8)
		# find all of the unique pert_ids in the data
		perts = [x.split(':::')[0].split('::')[1] for x in sco.pid] #perts is pert_desc
		pert_ids = [x.split(':')[1] for x in sco.pid]
		# unique_perts = set(perts)
		unique_perts = set(pert_ids)
		ctl_perts = []
		for i, unique_pert in enumerate(unique_perts):
			if unique_pert == 'DMSO':
				ctl_perts.append(unique_pert)
		unique_perts.difference_update(set(ctl_perts))
		#make pairing of pert id and pert_desc
		IDdescMatch = []
		for i, unique_pert in enumerate(unique_perts):
			IDdescMatch.append([unique_pert])
			i1 = pert_ids.index(unique_pert)
			desc1 = perts[i1]
			IDdescMatch[i].append(desc1)
		# grab the dose information
		dose = [float(x.split('::')[0].split(':')[2]) for x in sco.pid]
		# grab pert_descs
		desc = [x.split('::')[1].split(':::')[0] for x in sco.pid]
		# write sc plots to file
		num_perts = len(unique_perts)
		for i,unique_pert in enumerate(unique_perts):
			prog.update('making SC plots',i,num_perts)
			sco.plot(include=unique_pert,size=dose,title=(IDdescMatch[i][1] + ' - ' + unique_pert),pos_con=['None'],out=os.path.join(work_dir,'_'.join([unique_pert.replace(':','_'),'SC.png'])))
		### correlation matrix for each compound of interest
		# BRD-K01737880 - compound of interest 
		mPertDesc = ['AA09', 'AZD-1152HQPA', 'DOS', 'GSK-1070916', 'MLN8054', 'VX-680']
		mPertID = ['BRD-K29830875', 'AZD-1152HQPA', 'BRD-K01737880', 'GSK-1070916', 'MLN8054', 'VX-680']
		# icmpd1 = _all_indices(mPertID[1], db.get_column_meta('pert_id')) #grab the first compounds set of doses - are these the same for all compounds?
		value = mPertID[0]
		qlist = db.get_column_meta('pert_id')
		indices = []
		indx = -1
		while True:
			try:
				indx = qlist.index(value, indx+1)
				indices.append(indx)
			except ValueError:
				break
		uniqueDoses = [db.get_column_meta('pert_dose')[i] for i in indices]
		unqDoses = [float(x) for x in uniqueDoses]
		iSortDose = numpy.argsort(unqDoses) #order of doses sorted
		srtDoses = [uniqueDoses[i] for i in iSortDose]
		### make correlation matrix
		qPert = db.get_column_meta('pert_desc')
		qPertID = db.get_column_meta('pert_id')
		qDose = db.get_column_meta('pert_dose')
		qID = db.get_cids()
		corrMtrx = numpy.zeros((len(srtDoses),len(mPertID),len(mPertID))) #to become the matrix of pairwise correlations
		for id,dose1 in enumerate(srtDoses):
			# dose1 = srtDoses[0]
			indexDict = {} #for a single dose save the compound index
			for cmpdID in mPertID:
				for i in range(len(qPert)):
					if qDose[i] == dose1 and qPertID[i] == cmpdID:
						indexDict[cmpdID] = i
			# doseMtrx = numpy.zeros((db.matrix.shape[0],len(mPertID)))
			# for j,cmpdID in enumerate(mPertID):
			# 	ind1 = indexDict[cmpdID]
			# 	if not cmpdID == qPertID[ind1]:
			# 		print 'compound IDs do not match in corr matrix asignment'
			# 	col1 = db.matrix[:,ind1]
			# 	doseMtrx[:,j] = col1
			for i,cmd1 in enumerate(mPertID):
				for j,cmd2 in enumerate(mPertID):
					ind1 = indexDict[cmd1]
					ind2 = indexDict[cmd2]
					es = rslt.matrix[ind1,ind2] #run spearman correlation
					corrMtrx[id,i,j] = es
			plt.imshow(corrMtrx[id,:,:],interpolation='nearest',vmin=-1,vmax=1)
			plt.xticks( numpy.arange(6), mPertDesc, rotation=70)
			plt.yticks( numpy.arange(6), mPertDesc,)
			plt.title('pairwise signature correlation - ' + cell + ' ' +tim + ' dose=' + dose1)
			plt.colorbar()
			fname = cell + '_' +tim + '_dose_' + dose1 + '_pairwise_corr.png'
			outf = os.path.join(work_dir,fname)
			plt.savefig(outf,bbox_inches='tight', pad_inches=0)
			plt.close()
		# #show pairwise correlation of BRD-K01737880 with other compounds
		iSelect = 2 #pick the index for the compound of interest
		cmpdSelect = mPertID[iSelect]
		compareMtrx = numpy.zeros((len(srtDoses),len(mPertID)))
		for iCompare in range(len(mPertID)):
			doseLst = []
			for dose in range(len(srtDoses)):
				doseLst.append(corrMtrx[dose,iSelect,iCompare])
			plt.plot(doseLst) #plot each pairwise corr
			plt.ylim(ymin=-1.1)
			plt.ylim(ymax=1.1)
		plt.xlabel('dose')
		plt.ylabel('pairwise spearman corr')
		legLabel = []
		for x in mPertID:
			legLabel.append(x + ' - ' + cmpdSelect)
		pylab.legend((legLabel), shadow = True, loc = (0.01, 0.55))
		ltext = pylab.gca().get_legend().get_texts()
		pylab.setp(ltext[0], fontsize = 10)
		pylab.setp(ltext[1], fontsize = 10)
		pylab.setp(ltext[2], fontsize = 10)
		pylab.setp(ltext[3], fontsize = 10)
		fname = cell + '_' +tim + '_pairwise_corr_with_' + cmpdSelect +'.png'
		outf = os.path.join(work_dir,fname)
		plt.savefig(outf,bbox_inches=0, pad_inches=0)
		plt.close()

### run query internal to dataset - calculate pairwise ES score
#used: /xchip/cogs/hogstrom/scripts/PRISM/run_internal_PRISM_query_tool.sh
#load in results
cellLst = ['PC3', 'A375', 'MCF7']
timeLst = ['6H', '24H']
cell = 'MCF7'
tim = '24H'
cellLine = cell
timeP = tim
refControl = 'pc'
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/es/%s_%s_%s' % (cell,timeP,refControl)
if not os.path.exists(work_dir):
	os.mkdir(work_dir)
rsltIn = '/xchip/cogs/hogstrom/analysis/scratch/prism/MCF7/24H/wteslm/feb05/dose_analysis_tool.1360084956001/feb05/my_analysis.query_tool.2013020516300665/result_WTESLM.COMBINED_n59x59.gctx'
rslt = gct.GCT()
rslt.read(rsltIn)
#plot whole matrix
# mtrx = rslt.matrix
# plt.imshow(mtrx,interpolation='nearest')
# plt.colorbar()
# plt.show()

cids = rslt.get_cids()
# rPertIDs = [(x.split('_')[0].split(':')[2]) for x in sco.pid]
rPertIDs = [(x.split('_')[0]) for x in cids]
rDoses = [(x.split('_')[1].split('U')[0]) for x in cids]

mPertDesc = ['AA09', 'AZD-1152HQPA', 'DOS', 'GSK-1070916', 'MLN8054', 'VX-680']
mPertID = ['BRD-K29830875', 'AZD-1152HQPA', 'BRD-K01737880', 'GSK-1070916', 'MLN8054', 'VX-680']
# icmpd1 = _all_indices(mPertID[1], db.get_column_meta('pert_id')) #grab the first compounds set of doses - are these the same for all compounds?
value = mPertID[0]
qlist = rPertIDs
indices = []
indx = -1
while True:
	try:
		indx = qlist.index(value, indx+1)
		indices.append(indx)
	except ValueError:
		break
uniqueDoses = [rDoses[i] for i in indices]
unqDoses = [float(x) for x in uniqueDoses]
iSortDose = numpy.argsort(unqDoses) #order of doses sorted
srtDoses = [uniqueDoses[i] for i in iSortDose]
### make correlation matrix
qPert = rPertIDs
qPertID = rPertIDs
qDose = rDoses
corrMtrx = numpy.zeros((len(srtDoses),len(mPertID),len(mPertID))) #to become the matrix of pairwise correlations
for id,dose1 in enumerate(srtDoses):
	# dose1 = srtDoses[0]
	indexDict = {} #for a single dose save the compound index
	for cmpdID in mPertID:
		for i in range(len(qPert)):
			if qDose[i] == dose1 and qPertID[i] == cmpdID:
				indexDict[cmpdID] = i
	# doseMtrx = numpy.zeros((db.matrix.shape[0],len(mPertID)))
	for j,cmpdID in enumerate(mPertID):
		ind1 = indexDict[cmpdID]
		if not cmpdID == qPertID[ind1]:
			print 'compound IDs do not match in corr matrix asignment'
		col1 = db.matrix[:,ind1]
		doseMtrx[:,j] = col1
	for ind1 in range(len(mPertID)):
		for ind2 in range(len(mPertID)):
			col1 = doseMtrx[:,ind1]
			col2 = doseMtrx[:,ind2]
			r, p = scipy.stats.spearmanr(col1,col2) #run spearman correlation
			corrMtrx[id,ind1,ind2] = r

	plt.imshow(corrMtrx[id,:,:],interpolation='nearest',vmin=-1,vmax=1)
	plt.xticks( numpy.arange(6), mPertDesc, rotation=70)
	plt.yticks( numpy.arange(6), mPertDesc,)
	plt.title('pairwise signature correlation - ' + cell + ' ' +tim + ' dose=' + dose1)
	plt.colorbar()
	fname = cell + '_' +tim + '_dose_' + dose1 + '_pairwise_corr.png'
	outf = os.path.join(work_dir,fname)
	plt.savefig(outf,bbox_inches='tight', pad_inches=0)
	plt.close()
# #show pairwise correlation of BRD-K01737880 with other compounds
iSelect = 2 #pick the index for the compound of interest
cmpdSelect = mPertID[iSelect]
compareMtrx = numpy.zeros((len(srtDoses),len(mPertID)))
for iCompare in range(len(mPertID)):
	doseLst = []
	for dose in range(len(srtDoses)):
		doseLst.append(corrMtrx[dose,iSelect,iCompare])
	plt.plot(doseLst) #plot each pairwise corr
	plt.ylim(ymax=1.1)
plt.xlabel('dose')
plt.ylabel('pairwise spearman corr')
legLabel = []
for x in mPertID:
	legLabel.append(x + ' - ' + cmpdSelect)
pylab.legend((legLabel), shadow = True, loc = (0.01, 0.55))
ltext = pylab.gca().get_legend().get_texts()
pylab.setp(ltext[0], fontsize = 10)
pylab.setp(ltext[1], fontsize = 10)
pylab.setp(ltext[2], fontsize = 10)
pylab.setp(ltext[3], fontsize = 10)
fname = cell + '_' +tim + '_pairwise_corr_with_' + cmpdSelect +'.png'
outf = os.path.join(work_dir,fname)
plt.savefig(outf,bbox_inches=0, pad_inches=0)
plt.close()

# #make histogram of data
# corrLst = []
# for x in range(len(srtDoses)):
# 	for y in range(len(mPertID)):
# 		for z in range(len(mPertID)):
# 			corrLst.append(corrMtrx[x,y,z])
# plt.hist(corrLst,30)
# plt.show()
# #for each pairwise combination, make an x-y plot with correlations at dose
# for y in range(len(mPertID)):
# 	for z in range(len(mPertID)):
# 		corrByDose = []
# 		for x in range(len(srtDoses)):
# 			corrByDose.append(corrMtrx[x,y,z])
# 		plt.plot(corrByDose)
# plt.show()


### probe consistancy for each compound of interest
