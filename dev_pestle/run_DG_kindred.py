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

work_dir = '/xchip/cogs/hogstrom/analysis/informer_CTD/7May2013'

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
	#write drug signatures by cell line to a file
	outdir = os.path.join(work_dir,cell1)
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	sigF = os.path.join(outdir,cell1 + '_cp_sig_ids.grp')
	with open(sigF, 'w') as f:
		[f.write(x + '\n') for x in sigIDlist]
	# #write a gmt of the signatures
	# outdir = os.path.join(work_dir,cell1)
	# if not os.path.exists(outdir):
	# 	os.mkdir(outdir)
	# fup = os.path.join(work_dir,cell1,cell1 + '_sig_CPC_with_gene_target_up.gmt')
	# fdn = os.path.join(work_dir,cell1,cell1 + '_sig_CPC_with_gene_target_dn.gmt')
	# for sig in sigIDlist:
	# 	CM = mutil.CMapMongo()
	# 	mtchLst = CM.find({'sig_id':sig},{})
	# 	for match in mtchLst:
	# 		sig_id = str(match['sig_id'])
	# 		up50 = match['up50_lm']
	# 		dn50 = match['dn50_lm']
	# 		with open(fup,'a') as f:
	# 			f.write(sig_id + '\t' + sig_id + '\t')
	# 			for pt in up50:
	# 				f.write(pt + '\t')
	# 			f.write('\n')
	# 		with open(fdn,'a') as f:
	# 			f.write(sig_id + '\t' + sig_id + '\t')
	# 			for pt in dn50:
	# 				f.write(pt + '\t')
	# 			f.write('\n')

# generate the query command
for cell1 in uniqueLines:
	cellDir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	# sigF = os.path.join(cellDir, cell1 + '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
	cidF = glob.glob(cellDir + '/' + cell1 + '_all_CGS_sig_ids_n*.grp')[0]
	sigF = os.path.join(cellDir,cell1 + '_cp_sig_ids.grp')
	cmd = ' '.join(
			['rum -q local sig_query_tool',
			 '--sig_id ' + sigF,
			 '--metric wtcs',
			 '--column_space custom',
			 '--cid ' + cidF,
			 '--outdir ' + outdir,
			 '--mkdir false',
 			 '--mkdir save_tail',
			 '--score2rank_direc'])
	# os.system(cmd)



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

