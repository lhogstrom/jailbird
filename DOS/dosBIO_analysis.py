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

plate = 'DOSBIO001'
cellLine = 'PC3'
refControl ='pc'
timeP = '24H'
dataDir = '/xchip/cogs/projects/DOS/data'
gctfile = glob.glob(dataDir +'/%s/%s_%s_%s/by_pert_id_pert_dose/%s_%s_%s_COMPZ.MODZ_SCORE_LM_*x978.gctx' % (refControl,plate,cellLine,timeP,plate,cellLine,timeP))
gctfile = gctfile[0]

work_dir = '/xchip/cogs/projects/DOS/DOSBIO'
# work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/%s_%s_%s/fulcrum' % (cell,timeP,refControl)
# if not os.path.exists(work_dir):
# 	os.mkdir(work_dir)

gcto=gct.GCT()
gcto.read(gctfile)

sco = sc.SC()
sco.add_sc_from_gctx_meta(gctfile, verbose=False)
sco.set_thresh_by_specificity(0.8)
t1 = '%s_%s_%s' % (plate,cellLine,timeP)
outF = os.path.join(work_dir,'_'.join([t1,'SC.png']))
sco.plot(title=t1,out=outF)


pDescs = gcto.get_column_meta('pert_desc')
iBRDs = [i for i, x in enumerate(pDescs) if x[:3] == 'BRD']
ss= gcto.get_column_meta('distil_ss')
ssBRD = [ss[i] for i in iBRDs]
maxSS = max(ssBRD) 
iMaxSS = [i for i,x in enumerate(ss) if x == maxSS][0]
maxBRD = pDescs[iMaxSS]

brdCounts = []
fullPertList = []
for ipert in iBRDs:
	pert = pDescs[ipert]
	CM = mutil.CMapMongo()
	pert_List = CM.find({'pert_id':{'$regex':pert}},{'sig_id':True,'cell_id':True})
	if pert_List:
		brdCounts.append(len(pert_List))
		fullPertList.extend(pert_List)

cellsAll = [sig['cell_id'] for sig in fullPertList]
uniqCells = list(set(cellsAll))
### make gmt signature of drugs of interest
for cell1 in uniqCells:
	outdir = os.path.join(work_dir,cell1)
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	sigIDlist = []
	for sig in fullPertList:
		if sig['cell_id'] == cell1:
			sigIDlist.append(sig['sig_id'])
	sigIDlist = list(set(sigIDlist))
	#write drug signatures by cell line to a file
	sigF = os.path.join(outdir,cell1 + '_cp_sig_ids.grp')
	with open(sigF, 'w') as f:
		[f.write(x + '\n') for x in sigIDlist]
	#get all CGS for a cell line
	CM = mutil.CMapMongo()
	CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','cell_id':cell1},{'sig_id':True,'pert_iname':True})
	if CGSbyCell:
		nCGS = len(CGSbyCell)
		sigF = os.path.join(outdir, cell1+ '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
		with open(sigF, 'w') as f:
			for sig in CGSbyCell:
				f.write(sig['sig_id'] + '\n')

for cell1 in uniqCells:
	cellDir = os.path.join(work_dir,cell1) 
	cidF = glob.glob(cellDir + '/' + cell1 + '_all_CGS_sig_ids_n*.grp')
	if not cidF:
		continue
	cidF = cidF[0]
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	# sigF = os.path.join(cellDir, cell1 + '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
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
for icell, cell1 in enumerate(uniqCells):
	celldir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx'):
		print cell1 + 'no query result file'
		continue #if no results file, skip loop
	rsltF = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(rsltF)
	prog.update('analyzing {0}',icell,len(uniqCells))
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



