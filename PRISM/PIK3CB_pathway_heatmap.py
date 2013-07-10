 #! /usr/bin/env python
'''
examine/search for plates with dose data
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
import cmap.util.progress as progress
import subprocess
import datetime
import cmap.util.tool_ops as to
import cmap.analytics.dgo as dgo

work_dir = '/xchip/cogs/projects/DOS/10July_PI3K_pathway_heatmaps/KD_drug_wtcs'
metric='wtcs'
is_gold = False
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

# dpathwayList = ['PIK3CA', 'PIK3CD', 'PIK3CG', 'MTOR', 'AKT1', 'AKT2', 'PTEN']
pathwayList = ['PIK3C', 'MTOR', 'AKT', 'PTEN']
drug_list = ['BRD-K05756698','BRD-K12184916']

#cgs cell lines
CM = mu.CMapMongo()
CGSbyCell = CM.find({'pert_type':'trt_sh.cgs'},{'cell_id':True})
cgsCells = list(set(CGSbyCell))

#loop through and write cell IDs 
prog = progress.DeterminateProgressBar('genomic pert query')
for i,cell1 in enumerate(cgsCells):
    #get all CGS for a cell line
    prog.update('querying',i,len(cgsCells))
    CM = mu.CMapMongo()
    for gene in pathwayList:
        CGSbyCell = CM.find({'pert_iname':{'$regex':gene},'pert_type':'trt_sh.cgs','is_gold':True,'cell_id':cell1},{'sig_id':True,'pert_iname':True})
        if CGSbyCell:
            outdir = os.path.join(work_dir,cell1)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            nCGS = len(CGSbyCell)
            # sigF = os.path.join(outdir, cell1+ '_genomic_sig_ids_n' + str(nCGS) + '.grp')
            sigF = os.path.join(outdir, cell1+ '_genomic_sig_ids.grp')
            with open(sigF, 'a') as f:
                for sig in CGSbyCell:
                    f.write(sig['sig_id'] + '\n')
        else:
            print gene + ' in ' + cell1 + ' CGS not found'
    ### write drug sig ids
    for pert in drug_list: 
        CM = mu.CMapMongo()
        #grab one drug instance - dose > 5um
        if is_gold == True:
            pert_query = CM.find({'pert_id':pert,'is_gold':True,'cell_id':cell1},{'sig_id':True},limit=1)
        if is_gold == False:
            pert_query = CM.find({'pert_id':pert,'cell_id':cell1,'pert_dose':{'$gt':5}},{'sig_id':True},limit=1)
        if pert_query:
            with open(sigF, 'a') as f:
                for sig in pert_query:
                    f.write(sig+ '\n')
        else:
             print pert + ' in ' + cell1 + ' CGS not found'                  

### make query 
max_processes=10
processes = set()
for cell1 in cgsCells:
    cellDir = os.path.join(work_dir,cell1) 
    cidF = glob.glob(cellDir + '/' + cell1 + '_genomic_sig_ids.grp')
    if not cidF:
        continue
    cidF = cidF[0]
    outdir = os.path.join(work_dir,cell1,'sig_query_out')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # sigF = os.path.join(cellDir, cell1 + '_genomic_sig_ids_n' + str(nCGS) + '.grp')
    sigF = os.path.join(cellDir,cell1 + '_genomic_sig_ids.grp')
    cmd = ' '.join(['rum -q local -f sig_query_tool',
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
### dgo object 
# dg = dgo.QueryTargetAnalysis(out=work_dir)
# dg.make_result_frames(gp_type='KD',metric='wtcs')

### make result dataframe for each cell line - make heatmap
gp_type = 'KD'
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
    rsFrm = rslt.frame
    rsFrm = rsFrm.sort(axis=0)
    rsFrm = rsFrm.sort(axis=1)
    plt.pcolor(rsFrm,vmin=-1,vmax=1)
    igenes = [x.split(':')[1][:13] for x in rsFrm.index]
    cgenes = [x.split(':')[1][:13] for x in rsFrm.columns]
    plt.yticks(np.arange(0.5, len(rsFrm.index), 1), igenes)
    plt.xticks(np.arange(0.5, len(rsFrm.columns), 1), cgenes,rotation = 45)
    plt.title(cell1 + ' - CGS KD ' + metric)
    plt.colorbar()
    plt.savefig(os.path.join(work_dir,cell1+'_KD_heatmap.png'))
    plt.close()

