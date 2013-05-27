#! /usr/bin/env python
'''
analyze the DOSBIO plates - make sc plots etc
'''

import os
import cmap.io.gct as gct
import cmap.analytics.sc as sc
import glob as glob
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cmap.analytics.dose as doseClass
import subprocess
import os
import time

work_dir = '/xchip/cogs/projects/HOG/DG_connect'
#load in OMIM genes. Which ones have a CGS in > 4 cell lines? which ones are LM?
inFile = '/xchip/cogs/hogstrom/analysis/OMIM/OMIM_CGS.txt'
omimGeneList = []
with open(inFile,'rt') as f:
	for string in f:
		splt = string.split('\r')
		for i,line in enumerate(splt):
			if i == 0: # skip headder
				continue
			splt2 = line.split('\t')
			geneID = splt2[0] #the pert_id listed the line
			omimGeneList.append(geneID)

CM = mutil.CMapMongo()
CGSall = CM.find({'pert_type':'trt_sh.cgs'},{'sig_id':True,'pert_iname':True,'cell_id':True,'pert_id':True})
#which drugs to use --> informer set and HOG plate

### which genes have a CGS in > 4 cell lines
ominWithContext = []
for geneID in omimGeneList:
	cellLst = []
	sigIDLst = []
	for q in CGSall:
		if q['pert_iname'] == geneID:
			cellLst.append(q['cell_id'])
			sigIDLst.append(q['sig_id'])
	nCells = len(set(cellLst)) 
	if nCells > 4:
		ominWithContext.append(geneID)

### in which cell lines do we have more than 500 CGS sigs
cellFull = [q['cell_id'] for q in CGSall]
cgsCellsLong = list(set(cellFull))
cgsCells = []
for cell in cgsCellsLong:
	if cellFull.count(cell) > 500:
		cgsCells.append(cell)

### get HOG brds
file1 = '/xchip/obelix/pod/brew/pc/HOG001_A549_24H/by_pert_id_pert_dose/HOG001_A549_24H_COMPZ.MODZ_SCORE_LM_n288x978.gctx'
file2 = '/xchip/obelix/pod/brew/pc/HOG002_A549_24H/by_pert_id_pert_dose/HOG002_A549_24H_COMPZ.MODZ_SCORE_LM_n288x978.gctx'
hogPerts = []
#plate 1
dp = doseClass.DosePlate()
dp.add_from_gct(file1)
dp.examine_doses_tested()
hogPerts.extend(dp.perts_at_dose)
#plaste 2
dp = doseClass.DosePlate()
dp.add_from_gct(file2)
dp.examine_doses_tested()
hogPerts.extend(dp.perts_at_dose)
hogPertsShort = [x[:13] for x in hogPerts]

### pull sigs of hog perts from mongo
# for pert in hogPertsShort:
# for cell in cgsCells:
# 	CM = mutil.CMapMongo()
# 	cpSigs = CM.find({'pert_id':pert,'pert_type':'trt_cp'},{'sig_id':True,'cell_id':True,'pert_iname':True,'pert_id':True})
# 	sigLst = [x['sig_id'] for x in cpSigs]
# 	brdLst = [x['pert_id'] for x in cpSigs]
# 	cpSet = list(set(cpLst))
# 	cpSet.sort()
# 	brdSet = list(set(brdLst))

### make dictionary of pert_descs
pDescDict = {}
dp = doseClass.DosePlate()
dp.add_from_gct(file1)
for i,pert in enumerate(dp.pert_ids):
	if not pDescDict.has_key(pert[:13]):
		pDescDict[pert[:13]] = dp.pert_descs[i]
dp = doseClass.DosePlate()
dp.add_from_gct(file2)
for i,pert in enumerate(dp.pert_ids):
	if not pDescDict.has_key(pert[:13]):
		pDescDict[pert[:13]] = dp.pert_descs[i]

#plaste 2
dp = doseClass.DosePlate()
dp.add_from_gct(file2)

### organize sig ids by cell line
for cell1 in cgsCells:
	outdir = os.path.join(work_dir,cell1)
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	CM = mutil.CMapMongo()
	cpByCell = CM.find({'cell_id':cell1,'pert_type':'trt_cp'},{'sig_id':True,'cell_id':True,'pert_iname':True,'pert_id':True})
	sigIDlist = []
	pertIDlist = []
	for sig in cpByCell:
		if sig['pert_id'] in hogPertsShort:
			sigIDlist.append(sig['sig_id'])
			pertIDlist.append(sig['pert_id'])
	#count instances of pert in the cell line
	counts = dict()
	for i in pertIDlist:
	  counts[i] = counts.get(i, 0) + 1
	# cap at 10 signatures per compound per cell line
	iRemove = []
	for brd in counts:
		if counts[brd] > 10:
			icp = [i for i,x in enumerate(pertIDlist) if x == brd]
			irm = icp[11:] 
			iRemove.extend(irm)
	iKeep = set(range(len(sigIDlist)))
	iKeep = list(iKeep.difference(set(iRemove))) 
	sigListShort = [sigIDlist[i] for i in iKeep]
	### write pert and cgs sig IDs to file
	nSigs = len(sigListShort)
	sigF = os.path.join(outdir, cell1+ '_cp_sig_ids_n' + str(nSigs) + '.grp')
	with open(sigF, 'w') as f:
		[f.write(x + '\n') for x in sigListShort]
	#get all CGS for a cell line
	CM = mutil.CMapMongo()
	CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','cell_id':cell1},{'sig_id':True,'pert_iname':True})
	if CGSbyCell:
		nCGS = len(CGSbyCell)
		sigF = os.path.join(outdir, cell1+ '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
		with open(sigF, 'w') as f:
			for sig in CGSbyCell:
				f.write(sig['sig_id'] + '\n')

processes = set()
max_processes = 11
for cell1 in cgsCells:
	cellDir = os.path.join(work_dir,cell1) 
	cidF = glob.glob(cellDir + '/' + cell1 + '_all_CGS_sig_ids_n*.grp')
	if not cidF:
		continue
	cidF = cidF[0]
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	# sigF = os.path.join(cellDir, cell1 + '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
	sigF = glob.glob(cellDir + '/' + cell1 + '_cp_sig_ids_n*.grp')[0]
	cmd = ' '.join(['rum -q local sig_query_tool',
			 '--sig_id ' + sigF,
			 '--metric wtcs',
			 '--column_space custom',
			 '--cid ' + cidF,
			 '--out ' + outdir,
			 '--mkdir false',
			 '--save_tail false'])
	# os.system(cmd)
	processes.add(subprocess.Popen(cmd,shell=True))
	if len(processes) >= max_processes:
		os.wait()
		processes.difference_update(
			p for p in processes if p.poll() is not None)

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

#get unique brds with connection data
fullBRDs = []
for ind in df.index:
	brd = ind[0]
	fullBRDs.append(brd)
uniqBRDs = list(set(fullBRDs))
nCPs = []
for brd in uniqBRDs:
	cpRes = df.ix[brd]
	nCPs.append(cpRes.shape[0])
iMany = [i for i,x in enumerate(nCPs) if x > 12]
# for ipert in iMany:
for brd in uniqBRDs:
	# brd = uniqBRDs[ipert]
	cpRes = df.ix[brd]
	cpRank = dfRank.ix[brd]
	nCPs.append(cpRes.shape[0])
	meanSer = cpRes.mean()
	meanRnk = cpRank.mean()
	nullCnt = pd.isnull(cpRes)
	#how many cell lines were both the pert and target tested in
	valCounts = nullCnt.shape[0] - nullCnt.sum(axis=0)
	# CS_thresh = .45 #theshold for mean ss
	# for ind in meanSer[meanSer > CS_thresh].index:
	rank_thresh = 10 #theshold for mean ss
	for ind in meanRnk[meanRnk < rank_thresh].index:
		# print ind
		rnkSer = cpRank[ind]
		if min(rnkSer[rnkSer.notnull()]) < 1:
			if valCounts[ind] > 6:
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
				plt.savefig(os.path.join(work_dir,'cherry_pick',brd +'_' + ind + '_connections.png'))
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
				plt.savefig(os.path.join(work_dir,'cherry_pick',brd +'_' + ind + '_percent_rank.png'))
				plt.close()



###plot one instance in log form
for x in pDescDict:
	if pDescDict[x] == 'CD-1530':
		print x

for col in dfRank.columns:
	if col[:4] == 'MUC1':
		print col

brd = 'BRD-K25737009'
ind = 'MUC1_96H'
cpRes = df.ix[brd]
cpRank = dfRank.ix[brd]
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
		plt.scatter(np.log10(rnk),yVals)
plt.xlim((-2, 2))
plt.ylim((0,count+1))
plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
plt.xlabel('log10 (percent rank)')
plt.ylabel('cell line')
plt.title(brd + ' - ' + ind + ' connection - ' + gp_type)
plt.savefig(os.path.join(work_dir,'cherry_pick',brd +'_' + ind + '_log_rank.png'))
plt.close()

#raw percent rank
# brd = 'BRD-K25737009'
# ind = 'MUC1_96H'
brd = 'BRD-A35588707'
# ind = 'TOP2A_96H'
# ind = 'TOP2A_144H'
# ind = 'TOP2A_120H'
# ind = 'TOP1_96H'
ind = 'TOP1_144H'
cpRank = dfRank.ix[brd]
cpRes = df.ix[brd]
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
plt.savefig(os.path.join(work_dir,'cherry_pick',brd +'_' + ind + '_percent_rank1.png'))
plt.close()

### percent rank pvalues based on Navon et. al - Novel Rank-Based Statistical Methods Reveal MicroRNAs with Differential Expression in Multiple Cancer Types
# r- number of patients
rnkList = rnkList.dropna().values/100
r = len(rnkList)
tmpVec = []
for i,perRank in enumerate(rnkList):
	k=i+1
	tmpVec.append(np.power(perRank,k)*np.power(1-perRank,(r-k)))
# pval=sum(np.power(rnkList,range(1,r+1))*np.power(1-rnkList,r-range(1,r+1)))
pval=sum(tmpVec)

### check result
rFile = '/xchip/cogs/hogstrom/analysis/OMIM/cp_KD_connection/A549/sig_query_out/result_WTCS.LM.COMBINED_n1043x3620.gctx'
rslt = gct.GCT()
rslt.read(rFile)
cids = rslt.get_cids()
rids = rslt.get_rids()
#BRD-A15079084_CHERP_96H_connections
brd = 'BRD-A15079084'
kd = 'CHERP_96H'
ibrd = [i for i,x in enumerate(cids) if x.split(':')[1][:13] == brd]
ikd = [i for i,x in enumerate(rids) if x.split(':')[1][:5] == 'CHERP']
rslt.matrix[565,364]


## quick search for top2a
CM = mutil.CMapMongo()
geneFind = CM.find({'pert_type':'trt_sh.cgs','pert_iname':{'$regex':'TOP1'}},{'sig_id':True,'pert_iname':True,'cell_id':True,'pert_id':True})

# teniposide - 
#TOP2A - maybe positive connection
#TOP1 - maybe negative connection - PC3