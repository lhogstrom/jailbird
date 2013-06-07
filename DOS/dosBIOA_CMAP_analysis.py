#! /usr/bin/env python
'''
analyze the DOSBIO plates - combine with data in cmap database to perform query

use DOS signatures generate queries of the CGS data (cell line specific results)
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
import subprocess
import cmap.analytics.dgo as dgo

##make dosBIO gmts with matlab script
#/xchip/cogs/hogstrom/scripts/jailbird/DOS/make_dosBIO_gmt.m

##concatonate gmts for each cell line
work_dir = '/xchip/cogs/projects/DOS/6June2013b'
cellDirs = [f for f in os.listdir(work_dir)]
for cell in cellDirs:
	cellPath = work_dir + '/' + cell
	upFiles = glob.glob(cellPath + '/DOS*up50.gmt')
	dnFiles = glob.glob(cellPath + '/DOS*dn50.gmt')
	upFiles.sort()
	dnFiles.sort()
	outUP = os.path.join(cellPath,cell+'_up50.gmt')
	outDN = os.path.join(cellPath,cell+'_dn50.gmt')
	f = open(outUP, "a")
	for tempfile in upFiles:
		file1 = open(tempfile, "r")
		for line in file1:
			line = line.replace("_UP", "")
			f.write(line)
	f.close()
	f = open(outDN, "a")
	for tempfile in dnFiles:
		file1 = open(tempfile, "r")
		for line in file1:
			line = line.replace("_DN", "")
			f.write(line)
	f.close()

#make sc plots
cellLst = ['A375','A549','HT29','MCF7','PC3']
tpLst = ['24H']
plateLst = ['DOSBIO001','DOSBIO002']
refControl = 'pc'
dataDir = '/xchip/cogs/projects/DOS/data'
scDir = os.path.join(work_dir,'sc_plots')
if not os.path.exists(scDir):
    os.mkdir(scDir)
#make system call - run dose_plate_tool
for cellLine in cellLst:
	for timeP in tpLst:
		for plate in plateLst:
			gctfile = glob.glob(dataDir +'/%s/%s_%s_%s/by_pert_id_pert_dose/%s_%s_%s_COMPZ.MODZ_SCORE_LM_*x978.gctx' % (refControl,plate,cellLine,timeP,plate,cellLine,timeP))
			gctfile = gctfile[0]
			gcto=gct.GCT()
			gcto.read(gctfile)
			sco = sc.SC()
			sco.add_sc_from_gctx_meta(gctfile, verbose=False)
			sco.set_thresh_by_specificity(0.8)
			t1 = '%s_%s_%s' % (plate,cellLine,timeP)
			outF = os.path.join(scDir,'_'.join([t1,'SC.png']))
			sco.plot(title=t1,out=outF)

#get unique DOS compounds tested on plate
allBrds = [] #brds from both dosbio plates
for cellLine in cellLst:
	for timeP in tpLst:
		for plate in plateLst:
			gctfile = glob.glob(dataDir +'/%s/%s_%s_%s/by_pert_id_pert_dose/%s_%s_%s_COMPZ.MODZ_SCORE_LM_*x978.gctx' % (refControl,plate,cellLine,timeP,plate,cellLine,timeP))
			gctfile = gctfile[0]
			gcto=gct.GCT()
			gcto.read(gctfile)
			pDescs = gcto.get_column_meta('pert_desc')
			pID = gcto.get_column_meta('pert_id')
			iBRDs = [i for i, x in enumerate(pDescs) if x[:3] == 'BRD'] #only DOS compounds
			brdTested = [pID[i] for i in iBRDs]
			allBrds.extend(brdTested)
allBrds = list(set(allBrds))

#cell lines in which DOS were tested
CM = mutil.CMapMongo()
brdCells = CM.find({'sig_id':{'$regex':'DOS'}},{'cell_id':True})
dosCells = list(set(brdCells))
cellsTested = []
### get up/dn tags for the same cps from mongo - append to gmt
for cell in dosCells:
	#get all dose compounds recorded in that cell line
	CM = mutil.CMapMongo()
	pert_List = CM.find({'sig_id':{'$regex':'DOS'},'cell_id':cell},{'sig_id':True,'pert_id':True,'up50_lm':True,'dn50_lm':True})
	matchList = [x for x in pert_List if x['pert_id'] in allBrds]
	#skip cell line if not 
	if len(matchList) == 0:
		continue
	cellsTested.append(cell) #make new list of cells with both DOS cps of interest and CGSs tested
	cellPath = work_dir + '/' + cell
	if not os.path.exists(cellPath):
	    os.mkdir(cellPath)
	outUP = os.path.join(cellPath,cell+'_up50.gmt')
	outDN = os.path.join(cellPath,cell+'_dn50.gmt')
	for q in matchList:
		col_name = q['sig_id']
		#write to gmt list 
		with open(outUP,'a') as f:
			f.write(col_name + '\t' + col_name + '\t')
			for pt in q['up50_lm']:
				f.write(pt + '\t')
			f.write('\n')
		with open(outDN,'a') as f:
			f.write(col_name + '\t' + col_name + '\t')
			for pb in q['dn50_lm']:
				f.write(pb + '\t')
			f.write('\n')
	#get all CGS for a cell line
	CM = mutil.CMapMongo()
	CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','cell_id':cell},{'sig_id':True,'pert_iname':True})
	if CGSbyCell:
		nCGS = len(CGSbyCell)
		sigF = os.path.join(cellPath, cell+ '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
		with open(sigF, 'w') as f:
			for sig in CGSbyCell:
				f.write(sig['sig_id'] + '\n')

### run sig query
processes = set()
max_processes = 7
for cell1 in cellsTested:
    cellDir = os.path.join(work_dir,cell1) 
    tagUP = os.path.join(cellDir,cell1+'_up50.gmt')
    tagDN = os.path.join(cellDir,cell1+'_dn50.gmt')
    cidF = glob.glob(cellDir + '/' + cell1 + '_all_CGS_sig_ids_n*.grp')
    if not cidF:
        continue
    cidF = cidF[0]
    outdir = os.path.join(work_dir,cell1,'sig_query_out')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    cmd = ' '.join(['rum -q local sig_query_tool',
             '--uptag ' + tagUP,
             '--dntag ' + tagDN,
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


test1 = 'OEB001_A375_96H:BRDN0000399163:-666' #set random sig_id to initialize dgo object
test2 = 'OEB001_A375_96H:BRDN0000400484:-666'
reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir)
# dg.add_dictionary(targetDict=targetDict)
# dg.get_drug_kd_sig_ids()
# dg.run_drug_gene_query()
dg.make_result_frames()
dg.test_unknown_connections(make_graphs=True)

dg.FDR_correction(pDescDict=pDescDict)

### scratch 
#perform mongo query for each compound of interest
# brdCounts = []
# fullPertList = []
# for pert in allBrds:
# 	CM = mutil.CMapMongo()
# 	pert_List = CM.find({'pert_id':{'$regex':pert}},{'sig_id':True,'cell_id':True})
# 	if pert_List:
# 		brdCounts.append(len(pert_List))
# 		fullPertList.extend(pert_List)
# # where any of these cps tested on non 'DOS' plates
# plateNames = []
# for x in fullPertList:
# 	plateNames.append(x['sig_id'])
# notDosName = [x for x in plateNames if not x[:3] == 'DOS']

#check headers from up/dn gmt
# cellDirs = [f for f in os.listdir(work_dir)]
# for cell in cellDirs:
# 	cellPath = work_dir + '/' + cell
# 	outUP = os.path.join(cellPath,cell+'_up50.gmt')
# 	outDN = os.path.join(cellPath,cell+'_dn50.gmt')
# 	file1 = open(outUP, "r")
# 	headerListUP = []
# 	for line in file1:
# 		header = line.split('\t')[0]
# 		headerListUP.append(header)
# 	file1.close()
# 	file1 = open(outDN, "r")
# 	headerListDN = []
# 	for line in file1:
# 		header = line.split('\t')[0]
# 		headerListDN.append(header)
# 	file1.close()
# 	print cell
# 	print headerListUP == headerListDN


fig,ax = plt.subplots(1)
ax.plot(range(10))
ax.set_xscale('log')