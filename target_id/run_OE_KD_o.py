 #! /usr/bin/env python
'''
examining the connections between gene KD and OE
'''

import glob, HTML
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection as fdr
import cmap.util.mongo_utils as mu
from cmap.tools import sig_slice_tool
from cmap.io import gct,plategrp,rnk
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import subprocess
import datetime
import cmap.util.tool_ops as to
import random

metric = 'wtcs'

work_dir = '/xchip/cogs/projects/target_id/OE_KD_25June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

prog = progress.DeterminateProgressBar('perturbation cid query')
#cell lines in which OEs were recorded
CM = mu.CMapMongo()
allOE = CM.find({'pert_type':'trt_oe','is_gold':True},{'sig_id':True,'pert_iname':True,'cell_id':True})
cell_lines_tested = []
cellsAll = [sig['cell_id'] for sig in allOE]
uniqCells = list(set(cellsAll))

cell_lines_tested = []
BRDNdict = {} #dictionary of gene symbol for each BRDN number
for i,cell1 in enumerate(uniqCells):
    #get all CGS for a cell line
    prog.update('querying',i,len(uniqCells))
    CM = mu.CMapMongo()
    OEbyCell = CM.find({'pert_type':'trt_oe','is_gold':True,'cell_id':cell1},{'sig_id':True,'pert_iname':True})
    if OEbyCell:
        for query in OEbyCell:
            brdn = query['sig_id'].split(':')[1]
            BRDNdict[brdn] = query['pert_iname']
    if OEbyCell:
        CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','is_gold':True,'cell_id':cell1},{'sig_id':True,'pert_iname':True})
        cell_lines_tested.append(cell1)
        outdir = os.path.join(work_dir,cell1)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        nCGS = len(CGSbyCell)
        sigF = os.path.join(outdir, cell1+ '_KD_sig_ids_n' + str(nCGS) + '.grp')
        #write OE signatures by cell line to a file
        with open(sigF, 'w') as f:
            for sig in CGSbyCell:
                f.write(sig['sig_id'] + '\n')
        #write OE signatures by cell line to a file
        sigF = os.path.join(outdir,cell1 + '_OE_sig_ids.grp')
        with open(sigF, 'w') as f:
            for sig in OEbyCell:
                f.write(sig['sig_id'] + '\n')

### run query
max_processes = 10
processes = set()
for cell1 in uniqCells:
    cellDir = os.path.join(work_dir,cell1) 
    cidF = glob.glob(cellDir + '/' + cell1 + '_KD_sig_ids_n*.grp')
    if not cidF:
        continue
    cidF = cidF[0]
    outdir = os.path.join(work_dir,cell1,'sig_query_out')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # sigF = os.path.join(cellDir, cell1 + '_genomic_sig_ids_n' + str(nCGS) + '.grp')
    sigF = os.path.join(cellDir,cell1 + '_OE_sig_ids.grp')
    cmd = ' '.join(['rum -q local sig_query_tool',
             '--sig_id ' + sigF,
             '--metric ' + metric,
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

#which cell lines have a result dir
cellDirs = [f for f in os.listdir(work_dir) if os.path.isdir(work_dir+'/'+f)]
prog = progress.DeterminateProgressBar('dataframe read')
df = pd.DataFrame()
dfRank = pd.DataFrame()
#loop through each cell line add to df
for icell, cell1 in enumerate(cellDirs):
    #define directories and load in outputs
    outdir = os.path.join(work_dir,cell1,'sig_query_out')
    if not glob.glob(outdir + '/result_*.gctx'):
        print cell1 + ' no query result file'
        continue #if no results file, skip loop
    if metric == 'wtcs':
        rsltFile = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
    if metric == 'spearman':
        rsltFile = glob.glob(outdir + '/result_SPEARMAN_n*.gctx')[0]
    rslt = gct.GCT()
    rslt.read(rsltFile)
    prog.update(cell1,icell,len(cellDirs))
    rsltF = rslt.frame
    rsltF = rsltF.T
    indVals = rsltF.index.values
    brdnVals = [ind.split(':')[1] for ind in indVals]
    oeVals = [BRDNdict[brdn] for brdn in brdnVals]
    #make the column name gene and pert time
    geneVals = []
    for ind in rsltF.columns:
        # if gp_type == 'KD':
        gene = ind.split(':')[1]
        # if gp_type == 'OE':
        #     brdn = ind.split(':')[1]
        #     gene = self.BRDNdict[brdn]
        tp = ind.split(':')[0].split('_')[-1]
        gname = '_'.join([gene, tp])
        geneVals.append(gname)
    if len(geneVals) > len(set(geneVals)):
        print 'duplicate CGS for this celline'
    newF = rsltF
    newF.index = [oeVals, rsltF.index.values]
    # if gp_type == 'KD':
    newF.columns = geneVals
    # if gp_type == 'OE':
        # newF.columns = [geneVals, rsltF.columns.values]
    rankF = newF.rank(ascending=False,axis=1)
    perRankF = rankF / float(rankF.shape[1]) * 100.0
    #add cell line result to combined df
    if len(df) == 0:
        df = newF
        dfRank = perRankF
    else:
        df = pd.concat([df,newF],axis=0)
        dfRank = pd.concat([dfRank,perRankF],axis=0)
dfCS = df
dfRank = dfRank

## calculate gene-gene relationships
csValues = []
rnkValues = []
oeGenes = [ix[0] for ix in dfRank.index]
oeGeneset = set(oeGenes)
for gene in oeGeneset:
    kdTps = [x for x in dfRank.columns if x.split('_')[0] == gene]
    for kdTp in kdTps:
        csSer = dfCS.ix[gene][kdTp]
        csSer = csSer[csSer.notnull()]
        csValues.extend(csSer.values)
        rnkSer = dfRank.ix[gene][kdTp]
        rnkSer = rnkSer[rnkSer.notnull()]
        rnkValues.extend(rnkSer.values)
# plt.hist(rnkValues,30)
plt.subplot(211)
plt.hist(csValues,30)
plt.title('observed KD-OE - wtcs distribution')
# plt.hist(rnkValues,30)
# plt.title('observed KD-OE - percent rank distribution')

## calculate random gene-gene relationships
csRandom = []
rnkRandom = []
prog = progress.DeterminateProgressBar('perm')
for igene,gene in enumerate(oeGeneset):
    prog.update('gene',igene,len(oeGeneset))
    kdTps = [x for x in dfRank.columns if x.split('_')[0] == gene]
    if kdTps:
        for kdTp in kdTps:
            g = kdTp
            randKD = random.choice(dfRank.columns)
            csSer = dfCS.ix[gene][randKD]
            # print gene + ' - ' + randKD
            csSer = csSer[csSer.notnull()]
            csRandom.extend(csSer.values)
            rnkSer = dfRank.ix[gene][randKD]
            rnkSer = rnkSer[rnkSer.notnull()]
            rnkRandom.extend(rnkSer.values)
plt.subplot(212)
plt.title('Random KD-OE pairs - wtcs distribution')
plt.hist(csRandom,30)
plt.show()

#plot ranks
plt.subplot(211)
plt.hist(rnkValues,30)
plt.title('observed KD-OE - percent rank distribution')
plt.subplot(212)
plt.title('Random KD-OE pairs - percent rank distribution')
plt.hist(rnkRandom,30)
plt.show()

reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir)
# dg.add_dictionary(targetDict=targetDict)
dg.make_result_frames(gp_type='KD',metric='wtcs')

##
cellDirs = [f for f in os.listdir(work_dir) if os.path.isdir(work_dir+'/'+f)]
prog = progress.DeterminateProgressBar('dataframe read')
df = pd.DataFrame()
dfRank = pd.DataFrame()
#loop through each cell line add to df
for icell, cell1 in enumerate(cellDirs):
    #define directories and load in outputs
    outdir = os.path.join(work_dir,cell1,'sig_query_out')
    if not glob.glob(outdir + '/result_*.gctx'):
        print cell1 + ' no query result file'
        continue #if no results file, skip loop
    if metric == 'wtcs':
        rsltFile = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
    if metric == 'spearman':
        rsltFile = glob.glob(outdir + '/result_SPEARMAN_n*.gctx')[0]
    rslt = gct.GCT()
    rslt.read(rsltFile)
    prog.update(cell1,icell,len(cellDirs))
    rsltF = rslt.frame
    rsltF = rsltF.T
    indVals = rsltF.index.values
    pertVals = [ind.split(':')[1][:13] for ind in indVals]
    #make the column name gene and pert time
    geneVals = []
    for ind in rsltF.columns:
        if gp_type == 'KD':
            gene = ind.split(':')[1]
        if gp_type == 'OE':
            brdn = ind.split(':')[1]
            gene = self.BRDNdict[brdn]
        tp = ind.split(':')[0].split('_')[-1]
        gname = '_'.join([gene, tp])
        geneVals.append(gname)
    if len(geneVals) > len(set(geneVals)):
        print 'duplicate CGS for this celline'
    newF = rsltF
    newF.index = [pertVals, rsltF.index.values]
    if gp_type == 'KD':
        newF.columns = geneVals
    if gp_type == 'OE':
        newF.columns = [geneVals, rsltF.columns.values]
    rankF = newF.rank(ascending=False,axis=1)
    perRankF = rankF / float(rankF.shape[1]) * 100.0
    #add cell line result to combined df
    if len(df) == 0:
        df = newF
        dfRank = perRankF
    else:
        df = pd.concat([df,newF],axis=0)
        dfRank = pd.concat([dfRank,perRankF],axis=0)
self.dfCS = df
self.dfRank = dfRank