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
import HTML
import pandas as pd

work_dir = '/xchip/cogs/hogstrom/analysis/informer_CTD/9May2013'
if not os.path.exists(work_dir):
	os.mkdir(work_dir)

targetSheetF = '/xchip/cogs/hogstrom/analysis/informer_CTD/Informer_JDannot_v1.txt'
targetSheetXLS = '/xchip/cogs/hogstrom/analysis/informer_CTD/Informer_JDannot_v1.xlsx'
targetSheetF = '/xchip/cogs/hogstrom/analysis/informer_CTD/Informer_short.txt'
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
			if targets[0] == '' or targets[0] == '?':
				continue
			else:
				targetDict[pID] = targets
				pDescDict[pID] = pDesc

dgo = dgo.DGO()
dgo.add_dictionary(targetDict)

### make mongo query to establish cell lines with the apropriate pairs
#first check the cell lines which the hairpin are in
PlateLim = 'CPC'
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
		print 'no query result for ' + pert
		noCMAPmatch['compound'].append(pert)
	else:
		pertCellList = []
		pertSigList = []
		for sig in pert_List:
			pertCellList.append(sig['cell_id'])
			pertSigList.append(sig['sig_id'])
		for target in dgo.targetDict[pert]:
			CM = mutil.CMapMongo()
			#get cell lines tested for each CGS
			target_CGS_List = CM.find({'pert_iname':target,'pert_type':'trt_sh.cgs'},{'sig_id':True,'cell_id':True})
			if not target_CGS_List:
				print 'no query result for ' + target
				noCMAPmatch['gene'].append(target)
			else:
				cellList = []
				sigList = []
				for sig in target_CGS_List:
					cellList.append(sig['cell_id'])
					sigList.append(sig['sig_id'])
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

### write CGS sig_ids to grp file
for cell in uniqueLines:
	outdir = os.path.join(work_dir,cell)
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	CM = mutil.CMapMongo()
	CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','cell_id':cell},{'sig_id':True})
	nCGS = len(CGSbyCell)
	sigF = os.path.join(outdir, cell+ '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
	with open(sigF, 'w') as f:
		[f.write(x + '\n') for x in CGSbyCell]

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


### load in query result
for cell1 in uniqueLines:
	cellDir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	file_rslt = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	# file_rank = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(file_rslt)


prog = progress.DeterminateProgressBar('Drug-target')
for icell, cell1 in enumerate(cellList):
	celldir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]:
		print cell1 + 'no query result file'
		continue #if no results file, skip loop
	rsltF = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(rsltF)
	prog.update('analyzing {0}',icell,len(cell1))
	queryRids = rslt.get_rids()
	queryGenes = [x.split(':')[1] for x in queryRids]
	resCids = rslt.get_cids()
	resPerts = [x.split(':')[1] for x in resCids]
	#build rank matrix
	sortMatrix = np.argsort(rslt.matrix,axis=0)[::-1]
	rankMatrix = np.argsort(sortMatrix,axis=0) + 1
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
	targetRanks = {}
	targetPvals = {}
	targetCS = {}
	sumF = os.path.join(celldir,cell1 + '_drug-target_connection_summary.txt')
	headers = ['query_sig','target_KD_cgs','cs', 'query_rank', 'pval']
	with open(sumF,'w') as f:
		f.write('\t'.join(headers) + '\n')
		for pert in queryInd:
			targets = pertSigIDs[pert].keys()
			targetRanks[pert] = {}
			targetPvals[pert] = {}
			targetCS[pert] = {}
			for target in targets:
				iTarget = [i for i,x in enumerate(queryGenes) if x == target]
				if iTarget:
					tarRanks = []
					tarCS = []
					tarPvals = []
					for iTar in iTarget:
						iQuery = queryInd[pert]
						for iq in iQuery:
							cs = rslt.matrix[iTar,iq]
							rnk = rankMatrix[iTar,iq]
							pval = rnk/float(rslt.matrix.shape[0])
							tarCS.append(cs)
							tarRanks.append(rnk)
							tarPvals.append(pval)
							#write summary table
							sig = resCids[iq]
							f.write('\t'.join([sig,target,str(cs),str(rnk),str(pval)[:7]]) + '\n')
					targetCS[pert][target] = tarCS
					targetPvals[pert][target] = tarPvals
					targetRanks[pert][target] = tarRanks
	### put pvals into a vector, test for FDR
	pVec = []
	for pert in targetPvals:
		pvals = targetPvals[pert].values()
		for pval in pvals:
			pVec.extend(pval)
	#perform FDR on pVec
	[pBoolean, pCorrected] = fdr(pVec, alpha=0.1, method='indep')
	nPassFDR = sum(pBoolean)
	if nPassFDR:
		iPassFDR = [i for i,x in enumerate(pBoolean) if x == True]
		pPass = [pVec[i] for i in iPassFDR]
		pMaxThresh = max(pPass) #what is the largest pval that passed FDR
		#flag connections which pass FDR
		outF = os.path.join(celldir,cell1 + '_drug-target_FDR_summary.txt')
		headers = ['query_sig','target_KD_cgs','cs', 'query_rank', 'pval']
		with open(outF,'w') as f:
			f.write('\t'.join(headers) + '\n')
			for pert in targetPvals:
				qInds = queryInd[pert]
				for target in targetPvals[pert]:
					pvals = targetPvals[pert][target]
					for i,pval in enumerate(pvals):
						if pval <= pMaxThresh:
							query = resCids[qInds[i]]
							cs = targetCS[pert][target][i]
							rank = targetRanks[pert][target][i]
							#query, target gene, cs, rank, pval
							f.write('\t'.join([query,target,str(cs),str(rank),str(pval)[:7]]) + '\n')
	### check observed drug-target connections compared to null-distribution
	# of connections between drug and random gene knockdowns
	observedRanks = []
	observedCS = []
	for pert in targetRanks:
		for target in targetRanks[pert]:
			observedRanks.extend(targetRanks[pert][target][:]) #create list of observed ranks
			observedCS.extend(targetCS[pert][target][:]) #create list of observed ranks
	n_rand = 1000
	randMtrxRank = np.zeros((len(observedRanks),n_rand))
	randMtrxCS = np.zeros((len(observedRanks),n_rand))
	for perm in range(n_rand):
		randRnkList = []
		randCSList = []
		for pert in targetRanks:
			iQuery = queryInd[pert]
			for target in targetRanks[pert]:
				iRand = np.random.randint(rslt.matrix.shape[0], size=1)
				tarRanks = rankMatrix[iRand,iQuery]
				tarCS = rslt.matrix[iRand,iQuery]
				randRnkList.extend(tarRanks) 
				randCSList.extend(tarCS) 
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
			sortSumFr = sumFrame.sort(columns='cs',ascending=False)
			sortSumFr = sortSumFr.set_index([range(len(sortSumFr))])
			table_data = []
			for i in range(len(sortSumFr)):
				line = [sortSumFr.query_sig[i],
						sortSumFr.target_KD_cgs[i],
						str(sortSumFr.cs[i]),
						str(sortSumFr.query_rank[i]),
						str(sortSumFr.pval[i])]
				table_data.append(line[:])
			htmlcode = HTML.table(table_data)
			f.write(htmlcode + '\n')



### write connections not found in mongo



# check rank results
cids = rslt.get_cids()
rids = rslt.get_rids()
iARF = rids.index('CGS001_MCF7_96H:ARF1:2')
i743 = cids.index('CPC006_MCF7_24H:BRD-A31107743-001-01-3:0.09')
rslt.matrix[iARF,i743]
rankMatrix[iARF,i743]
### get ranks to work correctly

