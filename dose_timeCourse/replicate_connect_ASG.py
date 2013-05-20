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

work_dir = '/xchip/cogs/projects/ASG_dose_time/replicate_connect/scratch'
fname = '/xchip/obelix/pod/tnwork/ASG_combined/ASG_merged_ZSPCQNORM_n2841x978.gctx'
db = gct.GCT()
db.read(fname) 
nProbes = db.matrix.shape[0]


pert_descs = db.get_column_meta('pert_desc')
pertIDs = db.get_column_meta('pert_id')
cellIDs = db.get_column_meta('cell_id')
timePs = db.get_column_meta('pert_time')
doses = db.get_column_meta('pert_dose')
cids = db.get_cids()
unique_perts = list(set(pert_descs))

# ## identify replicates ###
# must match on the following criteria: cellID, timePoint, pertID, dose
# replicateDct = {}
# searchLst = range(len(pert_descs))
# #group the indices for each replicate - save to a dictionary
# while (len(searchLst) > 0):
# 	i = searchLst[0]
# 	for j in searchLst[1:]:
# 		if (cellIDs[i] == cellIDs[j]
# 				and timePs[i] == timePs[j]
# 				and pert_descs[i] == pert_descs[j] 
# 				and doses[i] == doses[j]):
# 			if replicateDct.has_key(cellIDs[i]+', '+timePs[i]+', '+pertIDs[i]+', '+doses[i]):
# 				#replicateDct[cellIDs[i]+', '+timePs[i]+', '+pertIDs[i]+', '+doses[i]].append(i)
# 				replicateDct[cellIDs[i]+', '+timePs[i]+', '+pertIDs[i]+', '+doses[i]].append(j)
# 			else:
# 				replicateDct[cellIDs[i]+', '+timePs[i]+', '+pertIDs[i]+', '+doses[i]] = [i]
# 				replicateDct[cellIDs[i]+', '+timePs[i]+', '+pertIDs[i]+', '+doses[i]] = [j] #create key in dictionary if it doesn't exist
# 			# remove these replicates from the search list
# 			searchLst.remove(j)
# 	searchLst.remove(i)		

#check replicate groups
# 'ASG001_MCF7_6H_X2_B7_DUO52HI53LO:N11' - has only three replicates?
# 'MCF7, 6.0, BRD-K89732114-300-06-3, 2.0'

# cell1 = 'MCF7'
# time1 = '6.0'
# cmpd1 = 'BRD-K89732114-300-06-3'
# dose1 = '2.0'
# repList1 = []
# for i in range(len(pert_descs)):
# 	if (cellIDs[i] == cell1
# 			and timePs[i] == time1
# 			and pertIDs[i] == cmpd1 
# 			and doses[i] == dose1):
# 		repList1.append(i)


replicateDct = {}
for l in set(cellIDs):
	for m in set(timePs):
		for n in set(pertIDs):
			for o in set(doses):
				cell1 = l
				time1 = m
				cmpd1 = n
				dose1 = o
				repList1 = []
				for i in range(len(pert_descs)):
					if (cellIDs[i] == cell1
							and timePs[i] == time1
							and pertIDs[i] == cmpd1 
							and doses[i] == dose1):
						repList1.append(i)
				if len(repList1) >= 1:
					replicateDct[l+', '+m+', '+n+', '+o] = repList1

#check distribution of replicates
LengthReps = []
for rep in replicateDct:
	lst1 = replicateDct[rep]
	LengthReps.append(len(lst1))

shrtLst = []
shrtLst = [x for x in LengthReps if x<20]
# import pylab
# def fig_1():
# 	n, bins, patches = pylab.hist(shrtLst, normed=1)
# 	l, = pylab.plot(bins, pylab.normpdf(bins, 0.0, 1.0), 'r--', label='fit', linewidth=3)
# 	pylab.legend([l, patches[0]], ['fit', 'hist'])




# # calcualte spearman corr among replicates
# # take one%pa replicate set
# for x in replicateDct:
# 	inds = replicateDct[x]
# 	#create z score matrix for replicates
# 	mtrx = numpy.zeros(shape=(nProbes,len(inds))
# 	for l,m in enumerate(inds):
# 		mtrx[:,l] = db.matrix[:,m]

#load in query_tool results
fname = '/xchip/cogs/projects/ASG_dose_time/replicate_connect/jan18/my_analysis.query_tool.2013011812161597/result_WTESLM.COMBINED_n2841x2841.gctx'
rslt = gct.GCT()
rslt.read(fname) 
rsltCids = rslt.get_cids()
rsltRids = rslt.get_rids()

### plots ###
nSmll = 100
smllMtx = rslt.matrix[0:nSmll,0:nSmll]
plt.imshow(smllMtx,interpolation='nearest',cmap='RdBu')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(rslt.matrix, interpolation='nearest')
fig.colorbar(cax)
plt.show()


### analyze replicate connections ### 
ESmat = rslt.matrix
iES = ESmat.argsort(axis=0)[::-1] #sort ascending
n_inst = len(iES[:,1])

repRnkAll = []
# repGrp = replicateDct['MCF7, 24.0, BRD-A60070924-001-01-8, 0.08']
#repGrp = replicateDct['PC3, 6.0, BRD-A19500257-001-04-7, 2.0']
repRnkDct = {}
#loop through all groups of replicates
for rep in replicateDct:
	repGrp = replicateDct[rep]
	grpID = [cids[i] for i in repGrp]
	GrpRnk = []
	for sID in grpID:
		d = rsltCids.index(sID)#find index of that instance in the results
		iE = iES[:,d] #ES sort index for one column
		sSigID = []
		for y in iE:
			sSigID.append(rsltRids[y]) #make sorted sig ID list
		cutGrp = grpID[:]
		cutGrp.remove(sID)
		repRnk = [sSigID.index(x) for x in cutGrp]
		GrpRnk.extend(repRnk)
		repRnkAll.extend(repRnk)
	repRnkDct[rep] = GrpRnk


dSet = list(set(doses))
dSetF = [float(x) for x in dSet]
dSetF.sort()

for l in set(cellIDs):
	for m in set(timePs):
		for n in set(pertIDs):
			# l = 'MCF7'
			# m = '6.0'
			# n = 'BRD-K89732114-300-06-3'
			ip = pertIDs.index(n)
			pDesc = pert_descs[ip] #find the pert_desc that matches pert_id
			fig = plt.figure()
			j=1
			# for o in set(doses):
			for o in dSetF:
				dose1 = str(o)
				if repRnkDct.has_key(l+', '+m+', '+n+', '+str(o)):
					repResult = repRnkDct[l+', '+m+', '+n+', '+str(o)]
					subPltN = 410 +j #set subplot location
					fig.add_subplot(subPltN)
					fig.subplots_adjust(hspace = 0.6)
					plt.hist(repResult,50,range =[0,len(iE)])
					plt.xlim(0, len(iE))
					plt.title('dose = ' + str(o)+'um')
					j = j + 1
			# plt.show()
			fname = l+'_'+m+'_'+pDesc+'_'+str(o)+'_replicate_connection.png'
			outf = os.path.join(work_dir,fname)
			print 'saving graph for ' + pDesc + '_' + dose1
			plt.savefig(outf, bbox_inches=0)

### write results to html file ###
if dSetF.__contains__(0.10):
	dSetF.remove(0.10)
for l in set(cellIDs):
	for m in set(timePs):
		topPercDic = {} #dict - key = pert_desc, value = top query percent for each dose
		for n in set(pertIDs):
			ip = pertIDs.index(n)
			pDesc = pert_descs[ip] #find the pert_desc that matches pert_id
			qTopDoses = []
			for i,o in enumerate(dSetF):
				dose1 = str(o)
				if repRnkDct.has_key(l+', '+m+', '+n+', '+str(o)):
					repResult = repRnkDct[l+', '+m+', '+n+', '+str(o)]
					nTop = len(iE)
					nPercent = .05 #percent of query list to look at
					cutOff = int(nTop*nPercent)
					TopCnt = [x for x in repResult if x < cutOff]
					nRepsInTop = len(TopCnt)/float(len(repResult))
					nRepsInTop = str(nRepsInTop)[:5] #convert to string and limit to 3 sig-figs
					qTopDoses.append(nRepsInTop)
			topPercDic[pDesc] = [qTopDoses]
		# add link to dose dictionary
		for x in set(pert_descs):
			if x == 'DMSO':
				continue
			else:
				lnk = l+'_'+m+'_'+x+'_'+str(o)+'_replicate_connection.png'
				topPercDic[x].append(lnk)
		if topPercDic.has_key('DMSO'):
			topPercDic.pop('DMSO')
		# setup jinga template 
		cmap_base_dir = '/'.join(os.path.dirname(cmap.__file__).split('/')[0:-1])
		env = jinja2.Environment(loader=jinja2.FileSystemLoader(cmap_base_dir + '/templates'))
		index_page_template = env.get_template('Link_color_table_Template.html')
		# index_links = [x + '_detail.html' for x in unique_pDescs]
		fout = l+'_'+m+'_'+'index.html'
		with open(os.path.join(work_dir,fout),'w') as f:
			f.write(index_page_template.render(title='Dose Analysis Results', data=topPercDic))




