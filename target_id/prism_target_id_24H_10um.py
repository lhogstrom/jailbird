 #! /usr/bin/env python

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo

import glob, HTML
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection as fdr
import cmap.util.mongo_utils as mu
from cmap.tools import sig_slice_tool
from cmap.io import gct,plategrp,rnk
import cmap.util.progress as progress
import subprocess

work_dir = '/xchip/cogs/projects/PRISM/target_id/19June2013_24H_10um'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

#load in gene/ cp list and make target dict
geneF = '/xchip/cogs/projects/PRISM/target_id/gene_targets.txt'
geneList = [line.strip() for line in open(geneF)]
## make target dict 
pDescDict = { "BRD-K59369769": "VX-680", "BRD-K36740062": "GSK-1070916","CMAP-AZD-1152": "AZD-1152HQPA","BRD-K29830875": "-666","BRD-K01737880": "-666","BRD-K83963101": "MLN8054"}
targetDict = {}
for pert in pDescDict:
	targetDict[pert] = geneList
# Query instances of cps
CM = mu.CMapMongo()
pert_List = CM.find({'pert_iname':{'$regex':'AZD-1152HQPA'}},{'pert_id':True,})
# erbb2Lst = CM.find({'pert_iname':{'$regex':'ERBB2'},'pert_type':'trt_sh.cgs'},{'sig_id':True,'pert_iname':True,'pert_id':True,})

### test KD
# reload(dgo)
dg = dgo.QueryTargetAnalysis(work_dir + '/drug_KD_connection')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD')





### do sig query that runs only 24H tp and 5 or 10um
genomic_pert='KD'
targetDict_loaded=True
pert_list=False
prog = progress.DeterminateProgressBar('perturbation cid query')
work_dir = dg.outputdir
if not os.path.exists(work_dir):
        os.mkdir(work_dir)
brdCounts = []
fullPertList = []
### for each drug perturbations of interest - find all instances in CMAP
if targetDict_loaded:
    for i,pert in enumerate(dg.targetDict):
        prog.update('querying cps',i,len(dg.targetDict))
        CM = mu.CMapMongo()
        # pert_query = CM.find({'pert_id':{'$regex':pert}},{'sig_id':True,'cell_id':True})
        #limit query to 5 or 10 um
        # pert_query = CM.find({'pert_id':{'$regex':pert},'pert_time':24,'pert_dose':10},{'sig_id':True,'cell_id':True})
        # pert_query = CM.find({'pert_id':{'$regex':pert},'pert_time':24,'pert_dose':{'$in': [5, 10]}},{'sig_id':True,'cell_id':True})
        pert_query = CM.find({'pert_id':{'$regex':pert},'pert_time':24,'pert_dose':{"$gt": 4.5}},{'sig_id':True,'cell_id':True})
        if pert_query:
            brdCounts.append(len(pert_query))
            fullPertList.extend(pert_query)
if pert_list:
    for i,pert in enumerate(pert_list):
        prog.update('querying cps',i,len(pert_list))
        CM = mu.CMapMongo()
        # pert_query = CM.find({'pert_id':{'$regex':pert}},{'sig_id':True,'cell_id':True})
        # pert_query = CM.find({'pert_id':{'$regex':pert},'pert_time':24,'pert_dose':{'$in': [5, 10]}},{'sig_id':True,'cell_id':True})
        pert_query = CM.find({'pert_id':{'$regex':pert},'pert_time':24,'pert_dose':{"$gt": 4.5}},{'sig_id':True,'cell_id':True})
        if pert_query:
            brdCounts.append(len(pert_query))
            fullPertList.extend(pert_query)
cell_lines_tested = []
cellsAll = [sig['cell_id'] for sig in fullPertList]
uniqCells = list(set(cellsAll))
prog = progress.DeterminateProgressBar('genomic pert query')
### 1) which celll lines have tested the target with a genomic pert - write the cids to a file
### 2) write cp cig ids to a file if there are CGSs in that cell line
if genomic_pert == 'OE':
    BRDNdict = {} #dictionary of gene symbol for each BRDN number
for i,cell1 in enumerate(uniqCells):
    #get all CGS for a cell line
    prog.update('querying',i,len(uniqCells))
    CM = mu.CMapMongo()
    if genomic_pert == 'KD':
        CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','cell_id':cell1},{'sig_id':True,'pert_iname':True})
    if genomic_pert == 'OE':
        CGSbyCell = CM.find({'pert_type':'trt_oe','cell_id':cell1},{'sig_id':True,'pert_iname':True})
        if CGSbyCell:
            for query in CGSbyCell:
                brdn = query['sig_id'].split(':')[1]
                BRDNdict[brdn] = query['pert_iname']
    if CGSbyCell:
        cell_lines_tested.append(cell1)
        outdir = os.path.join(work_dir,cell1)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        nCGS = len(CGSbyCell)
        sigF = os.path.join(outdir, cell1+ '_genomic_sig_ids_n' + str(nCGS) + '.grp')
        with open(sigF, 'w') as f:
            for sig in CGSbyCell:
                f.write(sig['sig_id'] + '\n')
        sigIDlist = []
        for sig in fullPertList:
            if sig['cell_id'] == cell1:
                sigIDlist.append(sig['sig_id'])
        sigIDlist = list(set(sigIDlist))
        #write drug signatures by cell line to a file
        sigF = os.path.join(outdir,cell1 + '_cp_sig_ids.grp')
        with open(sigF, 'w') as f:
            [f.write(x + '\n') for x in sigIDlist]
dg.cell_lines_tested = cell_lines_tested
if genomic_pert == 'OE':
    dg.BRDNdict = BRDNdict




dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
# dg.make_result_frames(gp_type='KD')
# dg.test_known_connections(gp_type='KD',pDescDict=pDescDict)
# # dg.test_unknown_rank_product(gp_type='KD')
# dg.FDR_correction(pDescDict=pDescDict)

dg.make_result_frames(gp_type='KD',metric='wtcs')
dg.test_known_connections(gp_type='KD',metric='wtcs',pDescDict=pDescDict)
dg.FDR_correction(pDescDict=pDescDict,metric='wtcs',outName='apriori_connections_pass_FDR',alpha=0.2)
dg.fdr_html_summary(fdrDir='apriori_connections_pass_FDR')
dg.gene_to_drug_similarity(testGene='PIK3CA',gp_type='KD',metric='wtcs',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)

#save variables for re-load
# dgCopy = dg
# dg.dfRank = dgCopy.dfRank
# dg.dfCS = dgCopy.dfCS
# dg.pVec = dgCopy.pVec
# dg.pDict = dgCopy.pDict

### TEST OE
# dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_OE_connection')
# dg.add_dictionary(targetDict=targetDictCGS)
# dg.get_sig_ids(genomic_pert='OE')
# # dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
# dg.make_result_frames(gp_type='OE')
# dg.test_known_OE_connections(pDescDict=pDescDict,gp_type='OE')
# dg.FDR_correction(pDescDict=pDescDict)

outpath = '/'.join([dg.outputdir,'apriori_connections_pass_FDR/'])
hyperLnkPath = '/xchip/cogs/web/icmap/hogstrom/PRISM_drug_gene_fdr'
cmd = ' '.join(['ln -s',
		 outpath,
		 hyperLnkPath])

