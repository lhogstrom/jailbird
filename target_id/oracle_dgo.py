import glob, HTML
import matplotlib.pyplot as plt
import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection as fdr
import cmap.util.mongo_utils as mu
from cmap.tools import sig_slice_tool
from cmap.io import gct,plategrp,rnk
import cmap.util.progress as progress
import subprocess
import cmap.util.tool_ops as to
import cmap.analytics.dgo as dgo

import unittest
import numpy as np
from cmap.analytics import oracle
from cmap.io import queryresult

work_dir = '/xchip/cogs/projects/target_id/oracle_16Jul'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)
#test query
pert_list = ['BRD-K69840642', 'BRD-K41859756']
pDescDict = {'BRD-K69840642':'ISOX','BRD-K41859756':'NVP-AUY922'}

dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.get_sig_ids(genomic_pert='KD',targetDict_loaded=False,pert_list=pert_list,is_gold=True)
dg.run_drug_gene_query(metric='spearman',max_processes=10)

### get query files
metric = 'spearman'
fileList = []
cellDirs = [f for f in os.listdir(dg.outputdir) if os.path.isdir(dg.outputdir+'/'+f)]
for icell, cell1 in enumerate(cellDirs):
    #define directories and load in outputs
    outdir = os.path.join(dg.outputdir,cell1,'sig_query_out')
    if not glob.glob(outdir + '/result_*.gctx'):
        print cell1 + ' no query result file'
        continue #if no results file, skip loop
    if metric == 'wtcs':
        rsltFile = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
    if metric == 'spearman':
        rsltFile = glob.glob(outdir + '/result_SPEARMAN_n*.gctx')[0]
    fileList.append(rsltFile)

# file = '/xchip/cogs/rogerhu/scratch/test_oracle/result_SPEARMAN_n14x10.gctx'
# file2 = '/xchip/cogs/projects/target_id/CTD2_25June2013/drug_KD_spearman/PC3/sig_query_out/result_SPEARMAN_n440x5163.gctx'
qres = queryresult.QueryResult()
qres.read_multiple(fileList,transpose=True,rank_axis=1)
ocl = oracle.Oracle(qres,out=dg.outputdir)
ocl.build_groups(["pert_iname"], ["pert_id"])
ocl.connect()
ocl.mk_summary_html()



files = fileList
rank_axis = 1
ds = pd.DataFrame()
for file in  files:
    dstmp = gct.fastGCT(file)
    rankTmp = dstmp.frame.rank(ascending=False,axis=rank_axis)
    pctrankTmp = rankTmp / float(rankTmp.shape[rank_axis])
    if len(ds) == 0:
        ds = dstmp.frame
        rank = rankTmp
        pctrank = pctrankTmp
    else:
        ds = pd.concat([ds,dstmp.frame],axis=0)
        rank = pd.concat([rank,rankTmp],axis=0)
        pctrank = pd.concat([pctrank,pctrankTmp],axis=0)
self.score = ds
self.rank = rank
self.pctrank = pctrank

### group

        if len(sample_fields)==1:
            sample_groupvar = self.query.rdesc.ix[:,sample_fields[0]].tolist()
        else:
            sample_groupvar = map(lambda x: '_'.join(x),
                                 self.query.rdesc.ix[:,sample_fields].values)
        self.row_groups = self.query.rdesc.index.groupby(np.array(sample_groupvar))
        
        if len(query_fields)==1:
            query_groupvar = self.query.cdesc.ix[:,query_fields[0]].tolist()
        else:
            query_groupvar = map(lambda x: '_'.join(x),
                                  query.cdesc.ix[:,query_fields].values)
        self.col_groups = self.query.cdesc.index.groupby(np.array(query_groupvar))


### connect 
qids='all'
sig_ids='all'
groupvar=None
method="rank_product",
irection='pos'
n_jobs=1
niter=10000
direction='pos'

ocl.method=method
assert direction in ['pos','neg','both'], "direction must be 'pos','neg',or 'both'"
if ocl.method not in CONNECTIVITY_FUNCTIONS:
    raise ValueError("Method '%s' not supported. " % ocl.method)
else:
    connector_class = CONNECTIVITY_FUNCTIONS[ocl.method]
    ocl.method_ = connector_class()
summary_list = []
for q in ocl.row_groups.keys():
    print "Connecting query group %s" % q
    for g in ocl.col_groups.keys():
        #print "Connecting to group %s"  % str(g)
        values = ocl.query.pctrank.ix[ocl.row_groups[q],ocl.col_groups[g]]
        entry = ocl.method_(ocl.query,values,q,g,direction)
        summary_list.append(entry)
        
ocl.connection_summary = pd.DataFrame(summary_list)
ocl.fdr_correction()
ocl.connection_summary.sort(columns='q_value',inplace=True)
ocl.connection_summary.to_csv(os.path.join(ocl.out,
                                            'connection_summary.txt'),
                               sep='\t',index=False)
print "All connections computed."



