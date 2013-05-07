#! /usr/bin/env python
'''
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

# work_dir = '/xchip/cogs/hogstrom/analysis/informer_CTD'
# #load each CGS file and find geneIDs that have a CGS
# CGSdir = '/xchip/cogs/projects/rnai_analysis/collapsed_signatures/Collapsed_Signatures_Landmark_978'
# CGSfiles = glob.glob(CGSdir + '/*_*_modz_by_pert_desc_signatures_any_target_n*x978.gctx')
# cellList = []
# for path1 in CGSfiles:
# 	name1 = path1.split('/')[-1]
# 	cell1 = name1.split('_')[0]
# 	cellList.append(cell1)
# geneIDdict = {} #save gene ids tested in each cell line
# for cell1 in cellList:
# 	gctFile = glob.glob(CGSdir + '/' + cell1 + '_*_modz_by_pert_desc_signatures_any_target_n*x978.gctx')
# 	cgs = gct.GCT()
# 	cgs.read(gctFile[0])
# 	geneIDdict[cell1] = cgs.get_column_meta('id')


# targetSheetF = '/xchip/cogs/hogstrom/analysis/informer_CTD/Informer_JDannot_v1.txt'
# targetSheetXLS = '/xchip/cogs/hogstrom/analysis/informer_CTD/Informer_JDannot_v1.xlsx'
# targetSheetF = '/xchip/cogs/hogstrom/analysis/informer_CTD/Informer_short.txt'
# targetDict = {}
# pDescDict = {}
# with open(targetSheetF,'rt') as f:
# 	for string in f:
# 		splt = string.split('\r')
# 		for i,line in enumerate(splt):
# 			splt2 = line.split('\t')
# 			pID = splt2[0] #the pert_id listed the line
# 			pDesc = splt2[1]
# 			targets = splt2[2]
# 			targets = targets.split(';')
# 			targetDict[pID] = targets
# 			pDescDict[pID] = pDesc

# # matchDict = {}
# # for pert in targetDict:
# # 	if targetDict[pert][0] == '' or targetDict[pert][0] == '?':
# # 		continue
# # 	else: #continue only if the drug has an expected target
# # 		print pert
# # 		print targetDict[pert]
# # 		#perform mongo query to see of there is 
# # 		CM = mutil.CMapMongo()
# # 		mtchLst = CM.find({'pert_id':{'$regex':pert}},{'pert_id':True,'cell_id':True})
# # 		for target in targetDict[pert]:
# # 			mtchLst = CM.find({'pert_desc':{'$regex':target}},{'pert_id':True,'cell_id':True,'pert_desc':True})

# # For a given gene target, in which cell lines do we have a colapsed gene sig?
# matchDict = {}
# for pert in targetDict:
# 	if targetDict[pert][0] == '' or targetDict[pert][0] == '?':
# 		continue
# 	else: #continue only if the drug has an expected target
# 		# print pert
# 		# print targetDict[pert]
# 		matchDict[pert] = {}
# 		for target in targetDict[pert]:
# 			#loop through each cell line and see if we have a cgs for the target
# 			matchDict[pert][target] = []
# 			for cell1 in geneIDdict.keys():
# 				isMatch = target in geneIDdict[cell1]
# 				if isMatch:
# 					matchDict[pert][target].append(cell1)
# #make a gmt gene sig file for each drug where there is a drug-target knockdown pair in a given cell line
# for cell1 in cellList:
# 	fup = os.path.join(work_dir,cell1,cell1 + '_sig_CPC006_target_paired_up.gmt')
# 	fdn = os.path.join(work_dir,cell1,cell1 + '_sig_CPC006_target_paired_dn.gmt')
# 	#pick one cell line - MCF7
# 	for pert in matchDict:
# 		print pert
# 		print matchDict[pert].keys()
# 		cellTested = False
# 		for gene in matchDict[pert].keys():
# 			if cell1 in matchDict[pert][gene]:
# 				cellTested = True
# 		if cellTested:
# 			# mtchLst = CM.find({'pert_id':{'$regex':pert},'cell_id':cell1},{'pert_id':True,'cell_id':True,'sig_id':True})
# 			CM = mutil.CMapMongo()
# 			mtchLst = CM.find({'pert_id':{'$regex':pert},'cell_id':cell1},{})
# 			if mtchLst:
# 				for match in mtchLst:
# 					sig_id = str(match['sig_id'])
# 					if sig_id[:3] == 'CPC': #continue only if a CPC plate
# 						up50 = match['up50_lm']
# 						dn50 = match['dn50_lm']
# 						with open(fup,'a') as f:
# 							f.write(sig_id + '\t' + sig_id + '\t')
# 							for pt in up50:
# 								f.write(pt + '\t')
# 							f.write('\n')
# 						with open(fdn,'a') as f:
# 							f.write(sig_id + '\t' + sig_id + '\t')
# 							for pt in dn50:
# 								f.write(pt + '\t')
# 							f.write('\n')
# #run sig_query_tool
# for cell1 in cellList:
# 	rankF1 = '/xchip/cogs/hogstrom/analysis/informer_CTD/' + cell1 + '/' + cell1 + '_modz_CGZ_RANK_n*x978.gctx'
# 	rankF = glob.glob(rankF1)[0]
# 	scoreF = glob.glob(os.path.join(work_dir,cell1,cell1+'_*_modz_by_pert_desc_signatures_any_target_n*x978.gctx'))[0]
# 	fup = os.path.join(work_dir,cell1,cell1 + '_sig_CPC006_target_paired_up.gmt')
# 	fdn = os.path.join(work_dir,cell1,cell1 + '_sig_CPC006_target_paired_dn.gmt')
# 	outdir = os.path.join(work_dir,cell1,'sig_query_out')
# 	if not os.path.exists(outdir):
# 		os.mkdir(outdir)
# 	cmd = 'rum -q local sig_query_tool --build_id custom --rank ' + rankF + ' --score ' + scoreF + ' --uptag ' + fup + ' --dntag ' + fdn + ' --metric wtcs --mkdir false --out ' + outdir
# 	os.system(cmd)

### load in results file, find gene-gene connections
# rsltF = work_dir + '/may01/my_analysis.sig_query_tool.2013050118284095/result_WTCS.LM.COMBINED_n103x3232.gctx' 
prog = progress.DeterminateProgressBar('Drug-target')
outdir = os.path.join(work_dir,'6May2013')
for icell, cell1 in enumerate(cellList):
	if not glob.glob(os.path.join(work_dir,cell1,'sig_query_out','result_WTCS.LM.COMBINED_n*.gctx')):
		continue #if no results file, skip loop
	rsltF = glob.glob(os.path.join(work_dir,cell1,'sig_query_out','result_WTCS.LM.COMBINED_n*.gctx'))[0]
	rslt = gct.GCT()
	rslt.read(rsltF)
	prog.update('analyzing {0}',icell,len(cell1))
	queryGenes = rslt.get_rids()
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
	for pert in matchDict:
		iQuery = []
		iQuery = [i for i,x in enumerate(resPerts) if x[:13] == pert]
		if iQuery:
			queryInd[pert] = iQuery
	targetRanks = {}
	targetPvals = {}
	targetCS = {}
	for pert in queryInd:
		targets = matchDict[pert].keys()
		targetRanks[pert] = {}
		targetPvals[pert] = {}
		targetCS[pert] = {}
		for target in targets:
			iTarget = []
			iTarget = [i for i,x in enumerate(queryGenes) if x == target]
			if iTarget:
				queryGenes.index(target) #find the indices of the target in the query result
				# iQuery = 
				iQuery = queryInd[pert]
				tarCS = rslt.matrix[iTarget,iQuery]
				tarRanks = rankMatrix[iTarget,iQuery]
				tarPvals = tarRanks/float(rslt.matrix.shape[0])
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
		outF = os.path.join(outdir,cell1 + '_drug-target_connection_summary.txt')
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
	plt.savefig(os.path.join(outdir,cell1 + '_drug_target_rank_dist.png'))
	plt.close()
	#plot CS scores
	flat = randMtrxCS.flatten()
	h1 = plt.hist(flat,50,color='b',range=[min(flat),max(flat)],label=['drug-random gene'],normed=True,alpha=0.5)
	h2 = plt.hist(observedCS,50,color='r',range=[min(flat),max(flat)],label='drug-target gene',normed=True,alpha=0.5)
	plt.legend()
	plt.xlabel('cs score')
	plt.ylabel('freq')
	plt.title(cell1 + ' - connection of drug and target gene knockdown')
	plt.savefig(os.path.join(outdir,cell1 + '_drug_target_CS_dist.png'))
	plt.close()

### make html summary page
indexF = os.path.join(outdir,'index.html')
with open(indexF,'w') as f:
	lineWrite = '<h2><CENTER>Connection between a drug and the knockdown of its gene target by cell line <CENTER></h2>'
	f.write(lineWrite + '\n')
	for cell1 in cellList:
		outF = os.path.join(outdir,cell1 + '_drug_target_CS_dist.png')
		if os.path.exists(outF):
			lineWrite =  '<a href="' + cell1 + '.html">' + cell1 + '</a> <BR>'
			f.write(lineWrite + '\n')
### make individual pages for each cell line
for cell1 in cellList:
	outF = os.path.join(outdir,cell1 + '_drug_target_CS_dist.png')
	if os.path.exists(outF):
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
			sigF = os.path.join(outdir,cell1 + '_drug-target_connection_summary.txt')
			if os.path.exists(sigF):
				table_data = []
				with open(sigF,'rt') as Fread:
					for line in Fread:
						table_data.append(line.split('\t'))
				htmlcode = HTML.table(table_data)
				f.write(htmlcode)
			else:
				f.write('<TD>No tested connections pass FDR</TD>')
			lineWrite = '<h2><CENTER>All connections tested <CENTER></h2>'
			f.write(lineWrite + '\n')
			### write connections tested


## what are the affogato connections strengths?


