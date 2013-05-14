#! /usr/bin/env python
'''
set files up to run 
check the connection between a drug and the knockdown of its gene target
'''

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mutil
import glob
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection as fdr
import cmap.util.progress as progress
import cmap.analytics.dgo as dgo
import cmap.io.rnk as rnk
import HTML
import pandas as pd

work_dir = '/xchip/cogs/hogstrom/analysis/informer_CTD/14May2013'
if not os.path.exists(work_dir):
	os.mkdir(work_dir)

targetSheetF = '/xchip/cogs/hogstrom/analysis/informer_CTD/Informer_JDannot_v1.txt'
targetSheetXLS = '/xchip/cogs/hogstrom/analysis/informer_CTD/Informer_JDannot_v1.xlsx'
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

dgo = dgo.DGO()
dgo.add_dictionary(targetDict)

# question, how to know which cell lines they pair with
#grab all genetic perturbations - find which cell lines they occured in
CM = mutil.CMapMongo()
GPlist = CM.find({'pert_type':'trt_sh.cgs'},{'sig_id':True,'cell_id':True})
cellList = [x['cell_id'] for x in GPlist]
uniqueLines = list(set(cellList))

### write sig_ids to grp file by cell line
### store sig ids by gene ID/ cell type 
sigGeneDict = {}
geneIDDict = {}
for cell in uniqueLines:
	sigGeneDict[cell] = {}
	geneIDDict[cell] = {}
	outdir = os.path.join(work_dir,cell)
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	CM = mutil.CMapMongo()
	CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','cell_id':cell},{'sig_id':True,'pert_iname':True})
	# CGSbyCell = CM.find({'pert_type':'trt_oe','cell_id':cell},{'sig_id':True})
	nCGS = len(CGSbyCell)
	sigF = os.path.join(outdir, cell+ '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
	with open(sigF, 'w') as f:
		for sig in CGSbyCell:
			f.write(sig['sig_id'] + '\n')
			sigGeneDict[cell][sig['sig_id']] = sig['pert_iname']
			if not geneIDDict.has_key(sig['pert_iname']):
				geneIDDict[sig['pert_iname']] = {}
			if geneIDDict[sig['pert_iname']].has_key(cell):
				geneIDDict[sig['pert_iname']][cell].append(sig['sig_id'])
			else:
				geneIDDict[sig['pert_iname']][cell] = [sig['sig_id']]

PlateLim = []
pertSigIDs = {}
CGSlines = {}
CGSsigID = {}
noCMAPmatch = {}
noCMAPmatch['compound'] = []
noCMAPmatch['gene'] = []
for pert in dgo.targetDict:
	pertSigIDs[pert] = {}
	CGSlines[pert] = {}
	CGSsigID[pert] = {}
	#get cell lines tested for each pert_id
	CM = mutil.CMapMongo()
	pert_List = CM.find({'pert_id':{'$regex':pert}},{'sig_id':True,'cell_id':True})
	if not pert_List:
		print 'no query result for ' + pert + ', ' + pDescDict[pert]
		noCMAPmatch['compound'].append(pert)
	else:
		pertCellList = []
		pertSigList = []
		for sig in pert_List:
			pertCellList.append(sig['cell_id'])
			pertSigList.append(sig['sig_id'])
		for target in dgo.targetDict[pert]:
			# CM = mutil.CMapMongo()
			# #get cell lines tested for each KD CGS
			# target_CGS_List = CM.find({'pert_iname':target,'pert_type':'trt_sh.cgs'},{'sig_id':True,'cell_id':True})
			#get cell lines tested for each OE sig
			# target_CGS_List = CM.find({'pert_iname':target,'pert_type':'trt_oe'},{'sig_id':True,'cell_id':True})
			if not geneIDDict.has_key(target):
				print 'no query result for ' + target
				noCMAPmatch['gene'].append(target)
			else:
				cellList = []
				sigList = []
				for sig in geneIDDict[target]:
					cellList = geneIDDict[target].keys()
					sigList = geneIDDict[target].values()
				cellSet = list(set(cellList))
				CGSlines[pert][target] = cellSet
				CGSsigID[pert][target] = sigList
				pertSigIDs[pert][target] = {}
				#check when a drug and target have both been tested in the same cell line
				for cell in cellSet:
					if cell in pertCellList:
						iperts = [i for i,x in enumerate(pertCellList) if x ==cell]
						if PlateLim:
							cellSigs = [pertSigList[i] for i in iperts if pertSigList[i][:3] == PlateLim]
						else:
							cellSigs = [pertSigList[i] for i in iperts]
						pertSigIDs[pert][target][cell] = cellSigs

### combine info on queries
cellLines = []
allSigIDsPert = []
for pert in pertSigIDs:
	for target in pertSigIDs[pert]:
		for cell in pertSigIDs[pert][target]:
			cellLines.append(cell)
			allSigIDsPert.extend(pertSigIDs[pert][target][cell])
uniqueLines = list(set(cellLines))
allSigIDsCGS = []
for pert in CGSsigID:
	for target in CGSsigID[pert]:
		allSigIDsCGS.extend(CGSsigID[pert][target])

### make gmt signature of drugs of interest
for cell1 in uniqueLines:
	sigIDlist = []
	for pert in pertSigIDs:
		for target in pertSigIDs[pert]:
			for cell2 in pertSigIDs[pert][target]:
				if cell2 == cell1:
					sigIDlist.extend(pertSigIDs[pert][target][cell2])
	sigIDlist = list(set(sigIDlist))
	#write drug signatures by cell line to a file
	outdir = os.path.join(work_dir,cell1)
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	sigF = os.path.join(outdir,cell1 + '_cp_sig_ids.grp')
	with open(sigF, 'w') as f:
		[f.write(x + '\n') for x in sigIDlist]

# generate the query command
for cell1 in uniqueLines:
	cellDir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	# sigF = os.path.join(cellDir, cell1 + '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
	cidF = glob.glob(cellDir + '/' + cell1 + '_all_CGS_sig_ids_n*.grp')[0]
	sigF = os.path.join(cellDir,cell1 + '_cp_sig_ids.grp')
	cmd = ' '.join(['rum -q local sig_query_tool',
			 '--sig_id ' + sigF,
			 '--metric wtcs',
			 '--column_space custom',
			 '--cid ' + cidF,
			 '--out ' + outdir,
			 '--mkdir false',
			 '--save_tail false'])
	os.system(cmd)

prog = progress.DeterminateProgressBar('Drug-target')
targetRanks = {}
targetRnkPercs = {}
targetCS = {}
for icell, cell1 in enumerate(cellList):
	celldir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx'):
		print cell1 + 'no query result file'
		continue #if no results file, skip loop
	rsltF = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(rsltF)
	prog.update('analyzing {0}',icell,len(cellList))
	queryRids = rslt.get_rids()
	# for KD:
	queryGenes = [x.split(':')[1] for x in queryRids]
	# for OE:
	# queryGenes = [x.split('_')[1] for x in queryRids]
	resCids = rslt.get_cids()
	resPerts = [x.split(':')[1] for x in resCids]
	#build rank matrix
	sortMatrix = np.argsort(rslt.matrix,axis=0)[::-1]
	rankMatrix = np.zeros_like(sortMatrix)
	for icol in range(rankMatrix.shape[1]):
		col = sortMatrix[:,icol]
		rankMatrix[col,icol] = np.arange(len(sortMatrix[:,icol])) + 1
	#pandas rank
	pRslt = rslt.frame
	rankFrame = pRslt.rank(ascending=False)
	### find the indices of each query
	queryInd = {}
	for pert in pertSigIDs:
		iQuery = []
		iQuery = [i for i,x in enumerate(resPerts) if x[:13] == pert]
		if iQuery:
			queryInd[pert] = iQuery
	targetRanks[cell1] = {}
	targetRnkPercs[cell1] = {}
	targetCS[cell1] = {}
	sumF = os.path.join(celldir,cell1 + '_drug-target_connection_summary.txt')
	headers = ['query_sig','cp_pert_desc','target_KD_cgs','cs', 'query_rank', 'RnkPerc']
	with open(sumF,'w') as f:
		f.write('\t'.join(headers) + '\n')
		for pert in queryInd:
			targets = pertSigIDs[pert].keys()
			targetRanks[cell1][pert] = {}
			targetRnkPercs[cell1][pert] = {}
			targetCS[cell1][pert] = {}
			for target in targets:
				iTarget = [i for i,x in enumerate(queryGenes) if x == target]
				if iTarget:
					tarRanks = []
					tarCS = []
					tarRnkPercs = []
					for iTar in iTarget:
						iQuery = queryInd[pert]
						for iq in iQuery:
							cs = rslt.matrix[iTar,iq]
							rnk = rankMatrix[iTar,iq]
							RnkPerc = rnk/float(rslt.matrix.shape[0])
							tarCS.append(cs)
							tarRanks.append(rnk)
							tarRnkPercs.append(RnkPerc)
							#write summary table
							sig = resCids[iq]
							f.write('\t'.join([sig,pDescDict[pert],target,str(cs),str(rnk),str(RnkPerc)[:7]]) + '\n')
					targetCS[cell1][pert][target] = tarCS
					targetRnkPercs[cell1][pert][target] = tarRnkPercs
					targetRanks[cell1][pert][target] = tarRanks
	### put RnkPercs into a vector, test for FDR
	pVec = []
	for pert in targetRnkPercs[cell1]:
		RnkPercs = targetRnkPercs[cell1][pert].values()
		for RnkPerc in RnkPercs:
			pVec.extend(RnkPerc)
	#perform FDR on pVec
	[pBoolean, pCorrected] = fdr(pVec, alpha=0.1, method='indep')
	nPassFDR = sum(pBoolean)
	if nPassFDR:
		iPassFDR = [i for i,x in enumerate(pBoolean) if x == True]
		pPass = [pVec[i] for i in iPassFDR]
		pMaxThresh = max(pPass) #what is the largest RnkPerc that passed FDR
		#flag connections which pass FDR
		outF = os.path.join(celldir,cell1 + '_drug-target_FDR_summary.txt')
		headers = ['query_sig','cp_pert_desc','target_KD_cgs','cs', 'query_rank', 'RnkPerc']
		with open(outF,'w') as f:
			f.write('\t'.join(headers) + '\n')
			for pert in targetRnkPercs[cell1]:
				qInds = queryInd[pert]
				for target in targetRnkPercs[cell1][pert]:
					RnkPercs = targetRnkPercs[cell1][pert][target]
					for i,RnkPerc in enumerate(RnkPercs):
						if RnkPerc <= pMaxThresh:
							query = resCids[qInds[i]]
							cs = targetCS[cell1][pert][target][i]
							rank = targetRanks[cell1][pert][target][i]
							#query, target gene, cs, rank, RnkPerc
							f.write('\t'.join([query,pDescDict[pert],target,str(cs),str(rank),str(RnkPerc)[:7]]) + '\n')
	### check observed drug-target connections compared to null-distribution
	# of connections between drug and random gene knockdowns
	observedRanks = []
	observedCS = []
	for pert in targetRanks[cell1]:
		for target in targetRanks[cell1][pert]:
			observedRanks.extend(targetRanks[cell1][pert][target][:]) #create list of observed ranks
			observedCS.extend(targetCS[cell1][pert][target][:]) #create list of observed ranks
	n_rand = 1000
	randMtrxRank = np.zeros((len(observedRanks),n_rand))
	randMtrxCS = np.zeros((len(observedRanks),n_rand))
	for perm in range(n_rand):
		randRnkList = []
		randCSList = []
		for pert in queryInd:
			for target in targetRanks[cell1][pert]:
				nTargs = queryGenes.count(target)
				iRand = np.random.randint(rslt.matrix.shape[0], size=nTargs)
				for iTar in iRand:
					iQuery = queryInd[pert]
					cs = rslt.matrix[iTar,iQuery]
					rnk = rankMatrix[iTar,iQuery]
					randCSList.extend(cs)
					randRnkList.extend(rnk)
		randMtrxRank[:,perm] = randRnkList
		randMtrxCS[:,perm] = randCSList
	#plot ranks
	flat = randMtrxRank.flatten()
	h1 = plt.hist(flat,30,color='b',range=[0,max(flat)],label=['drug-random gene'],normed=True,alpha=0.5)
	h2 = plt.hist(observedRanks,30,color='r',range=[0,max(flat)],label='drug-target gene',normed=True,alpha=0.5)
	plt.legend()
	plt.xlabel('query rank')
	plt.ylabel('freq')
	plt.title('connection of drug and target gene knockdown')
	plt.savefig(os.path.join(celldir,cell1 + '_drug_target_rank_dist.png'))
	plt.close()
	#plot CS scores
	flat = randMtrxCS.flatten()
	h1 = plt.hist(flat,50,color='b',range=[min(flat),max(flat)],label=['drug-random gene'],normed=True,alpha=0.5)
	h2 = plt.hist(observedCS,50,color='r',range=[min(flat),max(flat)],label='drug-target gene',normed=True,alpha=0.5)
	plt.legend()
	plt.xlabel('cs score')
	plt.ylabel('freq')
	plt.title(cell1 + ' - connection of drug and target gene knockdown')
	plt.savefig(os.path.join(celldir,cell1 + '_drug_target_CS_dist.png'))
	plt.close()

### drug-gene pairs across all cell lines
pairDict = {}
for pert in pertSigIDs:
	if pert in noCMAPmatch['compound']:
		continue
	else:
		for target in targetDict[pert]:
			if target in noCMAPmatch['gene']:
				continue
			else:
				key1 = pert + ':' + target
				pairDict[key1] = []
				for cell1 in targetRnkPercs:
					if pert in targetRnkPercs[cell1].keys():
						if target in targetRnkPercs[cell1][pert]:
							percRanks = targetRnkPercs[cell1][pert][target]
							pairDict[key1].extend(percRanks)
				### graph 1% club - WaddenGram
				if pairDict[key1]: 
					if min(pairDict[key1]) < .01 or max(pairDict[key1]) > .99:
						print key1
						sKeysStr = []
						for i,cell1 in enumerate(targetRnkPercs):
							sKeysStr.append(cell1)
							if pert in targetRnkPercs[cell1].keys():
								if target in targetRnkPercs[cell1][pert]:
									percRanks = targetRnkPercs[cell1][pert][target]
									yVals = np.repeat(i+1,len(percRanks))
									plt.scatter(percRanks,yVals)
						plt.xlim((0, 1))
						plt.ylim((0,i+2))
						plt.yticks(range(1, i + 2), sKeysStr, rotation = 0)
						plt.xlabel('percent rank of drug in CGS query')
						plt.ylabel('cell line')
						plt.title(pDescDict[pert] + ' - ' + target)
						plt.savefig(os.path.join(work_dir,pert +'_' + target + '_connections.png'))
						plt.close()

### make html summary page
indexF = os.path.join(work_dir,'index.html')
with open(indexF,'w') as f:
	lineWrite = '<h2><CENTER>Connection between a drug and the knockdown of its gene target by cell line <CENTER></h2>'
	f.write(lineWrite + '\n')
	for cell1 in cellList:
		cellDir = os.path.join(work_dir,cell1)
		if os.path.exists(cellDir):
			lineWrite =  '<a href="' + cell1 + '/' + cell1 + '.html">' + cell1 + '</a> <BR>'
			f.write(lineWrite + '\n')
### make individual pages for each cell line
for cell1 in cellList:
	outdir = os.path.join(work_dir,cell1)
	if os.path.exists(outdir):
		pageF = os.path.join(outdir,cell1+ '.html')
		with open(pageF,'w') as f:
			csF = os.path.join(outdir,cell1 + '_drug_target_CS_dist.png')
			lineWrite = '<BR>' + '<tr><td><img src=' + csF + '></td></tr>'
			f.write(lineWrite + '\n')
			rankF = os.path.join(outdir,cell1 + '_drug_target_rank_dist.png')
			lineWrite = '<BR>' + '<tr><td><img src=' + rankF + '></td></tr>'
			f.write(lineWrite + '\n')
			#generate and write summary table
			lineWrite = '<h2><CENTER>Connections which pass FDR (based on query rank) <CENTER></h2>'
			f.write(lineWrite + '\n')
			sigF = os.path.join(outdir,cell1 + '_drug-target_FDR_summary.txt')
			if os.path.exists(sigF):
				table_data = []
				with open(sigF,'rt') as Fread:
					for line in Fread:
						table_data.append(line.split('\t'))
				htmlcodeFDR = HTML.table(table_data)
				f.write(htmlcodeFDR + '\n')
			else:
				f.write('<TD>No tested connections pass FDR</TD>')
			lineWrite = '<h2><CENTER>All connections tested <CENTER></h2>'
			f.write(lineWrite + '\n')
			### write connections tested
			sumF = os.path.join(outdir,cell1 + '_drug-target_connection_summary.txt')
			sumFrame = pd.read_table(sumF)
			sortSumFr = sumFrame.sort(columns='RnkPerc',ascending=True)
			sortSumFr = sortSumFr.set_index([range(len(sortSumFr))])
			table_data = []
			for i in range(len(sortSumFr)):
				line = [sortSumFr.query_sig[i],
						sortSumFr.cp_pert_desc[i],
						sortSumFr.target_KD_cgs[i],
						str(sortSumFr.cs[i]),
						str(sortSumFr.query_rank[i]),
						str(sortSumFr.RnkPerc[i])]
				table_data.append(line[:])
			htmlcode = HTML.table(table_data)
			f.write(htmlcode + '\n')
			### write connections not found in mongo
			lineWrite = '<h2><CENTER>Perturbations not found in CMAP database <CENTER></h2>'
			f.write(lineWrite + '\n')
			tblNoFind = []
			for pert in noCMAPmatch['compound']:
				tblNoFind.append([pert])
			lineWrite = '<h3>Compounds </h3>'
			f.write(lineWrite + '\n')
			htmlcode = HTML.table(tblNoFind)
			f.write(htmlcode + '\n')
			lineWrite = '<h3>collapsed gene signatures </h3>'
			f.write(lineWrite + '\n')
			tblNoFind = []
			for pert in noCMAPmatch['gene']:
				tblNoFind.append([pert])
			htmlcode = HTML.table(tblNoFind)
			f.write(htmlcode + '\n')


### scratch ### 

# # check rank results
# cids = rslt.get_cids()
# rids = rslt.get_rids()
# iquery = queryGenes.index('RARB')
# hpID = queryRids[iquery]
# iARF = rids.index(hpID)
# i743 = cids.index('CPC006_VCAP_24H:BRD-K93176058-001-02-8:10')
# rslt.matrix[iARF,i743]
# rankMatrix[iARF,i743]
# col1  = rslt.matrix[:,i743]
# ### get ranks to work correctly
# array = np.array([4,2,7,1])
# temp = array.argsort()[::-1]
# ranks = np.empty(len(array), int)
# ranks[temp] = np.arange(len(array)) + 1

# sortMatrix = np.argsort(rslt.matrix,axis=0)[::-1]
# rankMatrix = np.zeros_like(sortMatrix)
# for icol in range(rankMatrix.shape[1]):
# 	col = sortMatrix[:,icol]
# 	rankMatrix[col,icol] = np.arange(len(sortMatrix[:,icol])) + 1

# # rankMatrix = np.argsort(sortMatrix,axis=0) + 1
# rnkRslt = rnk.rnk_objects_from_gctx(str(rsltF))
# rankMatrix = rnkRslt.matrix
