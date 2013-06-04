#! /usr/bin/env python

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mutil
import glob
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection as fdr
import cmap.util.progress as progress
import cmap.io.rnk as rnk
import HTML
import pandas as pd

work_dir = '/xchip/cogs/hogstrom/analysis/informer_CTD/14May2013'
### load informer drug-target data
targetSheetF = '/xchip/cogs/hogstrom/analysis/informer_CTD/Informer2_short.txt'
targetDict = {}
pDescDict = {}
with open(targetSheetF,'rt') as f:
	for string in f:
		splt = string.split('\r')
		for i,line in enumerate(splt):
			splt2 = line.split('\t')
			pID = splt2[0] #the pert_id listed the line
			pDesc = splt2[1]
			targets = splt2[2]
			targets = targets.split(';')
			if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
				continue
			else:
				targetDict[pID] = targets
				pDescDict[pID] = pDesc

# Create a pandas dataframe that lets you see connection results across 
# cell lines it is structured as follows:
# 	index1 = BRD short
# 	index2 = perurbation sig_id
# 	each column - a unique gene ID/ time point - representing the CGS for that gene, matching cell line
# 	cell line listed as a column
gp_type = 'KD' # genetic perturbation type
#which cell lines have a result dir
cellDirs = [f for f in os.listdir(work_dir) if os.path.isdir(work_dir+'/'+f)]
prog = progress.DeterminateProgressBar('Drug-target')
df = pd.DataFrame()
dfRank = pd.DataFrame()
#loop through each cell line add to df
# for icell, cell1 in enumerate(cgsCells):
for icell, cell1 in enumerate(cellDirs):
	#define directories and load in outputs
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx'):
		print cell1 + 'no query result file'
		continue #if no results file, skip loop
	rsltFile = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(rsltFile)
	prog.update('analyzing {0}',icell,len(cellDirs))
	rsltF = rslt.frame
	rsltF = rsltF.T
	indVals = rsltF.index.values
	pertVals = [ind.split(':')[1][:13] for ind in indVals]
	#make the column name gene and pert time
	geneVals = []
	for ind in rsltF.columns:
		gene = ind.split(':')[1]
		tp = ind.split(':')[0].split('_')[-1]
		gname = '_'.join([gene, tp])
		geneVals.append(gname)
	if len(geneVals) > len(set(geneVals)):
		print 'duplicate CGS for this celline'
	newF = rsltF
	newF.index = [pertVals, rsltF.index.values]
	newF.columns = geneVals
	rankF = newF.rank(ascending=False,axis=1)
	perRankF = rankF / float(rankF.shape[1]) * 100.0
	# newF['cell'] = cell1
	# perRankF['pert'] = pertVals
	# perRankF['cell'] = cell1
	#add cell line result to combined df
	if len(df) == 0:
		df = newF
		dfRank = perRankF
	else:
		df = pd.concat([df,newF],axis=0)
		dfRank = pd.concat([dfRank,perRankF],axis=0)

fullBRDs = []
for ind in df.index:
	brd = ind[0]
	fullBRDs.append(brd)
uniqBRDs = list(set(fullBRDs))
nCPs = []
for brd in uniqBRDs:
	cpRes = df.ix[brd]
	nCPs.append(cpRes.shape[0])
#loop through each drug-target pair
cols = dfRank.columns
pDict = {}
pVec = []
for brd in targetDict:
	#skip if not in query
	if not brd in uniqBRDs:
		continue
	targets = targetDict[brd]
	nTargets = []
	cpRes = df.ix[brd]
	cpRank = dfRank.ix[brd]
	nCPs.append(cpRes.shape[0])
	meanSer = cpRes.mean()
	meanRnk = cpRank.mean()
	nullCnt = pd.isnull(cpRes)
	#how many cell lines were both the pert and target tested in
	valCounts = nullCnt.shape[0] - nullCnt.sum(axis=0)
	for target in targets:
		tarList = [inst for inst in cols if inst.split('_')[0] == target]
		if len(tarList) == 0: #skip if drug target not tested
			continue
		for ind in tarList:
			print brd + ' - ' + ind
			rnkSer = cpRank[ind]/100
			rnkSer = rnkSer[rnkSer.notnull()]
			if len(rnkSer) == 0:
				continue
			### calculate p-value
			testStat = rnkSer.prod()
			n_obs = rnkSer.shape[0]
			#theoretical null
			n_rand = 10000
			prodDist = []
			for n in range(n_rand):
				percRnk = np.random.rand(1,n_obs)
				prodRnk = percRnk.prod()
				prodDist.append(prodRnk)
			#calculate p-values
			nullDist = np.array(prodDist)
			#number of null values more extreme than observed (one sided)
			exVals = nullDist[nullDist<testStat]
			nExtreme = len(exVals)
			pVal = (nExtreme+1)/float(len(nullDist))
			pVec.append(pVal)
			pDict[brd + '-' + ind] = pVal
			# plt.hist(prodDist,bins=np.logspace(-18, 0, 50))
			# plt.gca().set_xscale('log')
			# plt.show()
			# if valCounts[ind] > 9:
			if pVal < .01:
				print brd + ' ' + ind
				#cs wadden gram
				sKeysStr = []
				count = 0
				for i,cs in enumerate(cpRes[ind]):
					if pd.isnull(cs):
						continue
					else:
						count = count + 1
						sKeysStr.append(cpRes.index[i].split('_')[1])
						yVals = count
						plt.scatter(cs,yVals)
				plt.xlim((-1, 1))
				plt.ylim((0,count+1))
				plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
				plt.xlabel('wtcs')
				plt.ylabel('cell line')
				plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type)
				plt.savefig(os.path.join(work_dir,'drug_target_graphs',brd +'_' + ind + '_connections.png'))
				plt.close()
				#rank wadden gram
				sKeysStr = []
				count = 0
				rnkList = cpRank[ind]
				for i,rnk in enumerate(cpRank[ind]):
					if pd.isnull(rnk):
						continue
					else:
						count = count + 1
						sKeysStr.append(cpRes.index[i].split('_')[1])
						yVals = count
						plt.scatter(rnk,yVals)
				plt.xlim((0, 100))
				plt.ylim((0,count+1))
				plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
				plt.xlabel('percent rank')
				plt.ylabel('cell line')
				plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type)
				plt.savefig(os.path.join(work_dir,'drug_target_graphs',brd +'_' + ind + '_percent_rank.png'))
				plt.close()
### perform FDR correction
FDRtest = fdr(pVec)
pArr= np.array(pVec)
passed = pArr[FDRtest[0]]
pThreshFDR = max(passed)
cpPass = []
kdPass = []
print 'connections passing FDR correction:'
for test in pDict:
	if pDict[test] < pThreshFDR:
		print test
		drug = test.split('-')[0] + '-' + test.split('-')[1]
		cpPass.append(drug)
		kdPass.append(test.split('-')[-1])
#make graphs of connections that pass FRD 
for ibrd,brd in enumerate(cpPass):
	#skip if not in query
	if not brd in uniqBRDs:
		continue
	targets = targetDict[brd]
	cpRes = df.ix[brd]
	cpRank = dfRank.ix[brd]
	ind = kdPass[ibrd]
	#cs wadden gram
	sKeysStr = []
	count = 0
	for i,cs in enumerate(cpRes[ind]):
		if pd.isnull(cs):
			continue
		else:
			count = count + 1
			sKeysStr.append(cpRes.index[i].split('_')[1])
			yVals = count
			plt.scatter(cs,yVals)
	plt.xlim((-1, 1))
	plt.ylim((0,count+1))
	plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
	plt.xlabel('wtcs')
	plt.ylabel('cell line')
	plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type)
	plt.savefig(os.path.join(work_dir,'FDR_pass_graphs',brd +'_' + ind + '_connections.png'))
	plt.close()
	#rank wadden gram
	sKeysStr = []
	count = 0
	rnkList = cpRank[ind]
	for i,rnk in enumerate(cpRank[ind]):
		if pd.isnull(rnk):
			continue
		else:
			count = count + 1
			sKeysStr.append(cpRes.index[i].split('_')[1])
			yVals = count
			plt.scatter(rnk,yVals)
	plt.xlim((0, 100))
	plt.ylim((0,count+1))
	plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
	plt.xlabel('percent rank')
	plt.ylabel('cell line')
	plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type)
	plt.savefig(os.path.join(work_dir,'FDR_pass_graphs',brd +'_' + ind + '_percent_rank.png'))
	plt.close()

#calculate percent rank of individual drug-gene connections

#take product of percent ranks
#make histogram of this distribution
#calculate p-values

# based on random values
#simulate across
testRnks = np.random.rand(1,12)
testStat = testRnks.prod()
n_obs = testRnks.shape[1]
#theoretical null
n_rand = 10000
prodDist = []
for n in range(n_rand):
	percRnk = np.random.rand(1,n_obs)
	prodRnk = percRnk.prod()
	prodDist.append(prodRnk)
plt.hist(prodDist,bins=np.logspace(-18, 0, 50))
plt.gca().set_xscale('log')
plt.show()
#calculate p-values
nullDist = np.array(prodDist)
#number of null values more extreme than observed (one sided)
exVals = nullDist[nullDist<testStat]
nExtreme = len(exVals)
pVal = nExtreme/float(len(nullDist))


#how to make this a two tailed test --> want to see both positive/ negative

# what is a good way to reward results that are consistent within a cell line
# punish results that are inconsistent in a cell line
