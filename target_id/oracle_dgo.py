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
import cmap.analytics.dgo_oracle as dgo_oracle

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
targetDict = {'BRD-K69840642':['HSPA9'],'BRD-K41859756':['HSPA5']}

dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman2')
# dg.get_sig_ids(genomic_pert='KD',targetDict_loaded=False,pert_list=pert_list,is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)

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
reload(oracle)
ocl = oracle.Oracle(qres,out=dg.outputdir)
ocl.build_groups(["pert_iname","pert_time"], ["pert_id"])
ocl.connect_row()
ocl.connect_expected(targetDict)
ocl.connect(direction='both')

ocl.mk_summary_html()


### run oracle on CTD2 
work_dir = '/xchip/cogs/projects/target_id/CTD2_25June2013/oracle_17Jul'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_drug_targets.txt'
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
            targets = [x for x in targets if x != '']
            if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
                continue
            else:
                targetDict[pID] = targets
                pDescDict[pID] = pDesc

### scratch code playing with oracle modules ###


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


### connect - one row at a time
qids='all'
sig_ids='all'
groupvar=None
method="rank_product",
irection='pos'
n_jobs=1
niter=10000
direction='pos'

# assert self.row_groups is not None, "Must define row_groups. Call build_groups first."
# assert self.col_groups is not None, "Must define col_groups. Call build_groups first."
# assert direction in ['pos','neg','both'], "direction must be 'pos','neg',or 'both'"

# self.method=method
# self.method_ = method
summary_list = []
for q in ocl.row_groups.keys():
    print "Connecting query group %s" % q
    # row_vals = ocl.query.pctrank.ix[ocl.row_groups[q]]
    testStats = []
    counts = []
    prog = progress.DeterminateProgressBar('row test statistic')
    for ig,g in enumerate(ocl.col_groups.keys()):
        prog.update('querying cps',ig,len(ocl.col_groups.keys()))
        rnk_vals = ocl.query.pctrank.ix[ocl.row_groups[q],ocl.col_groups[g]]
        testStat = rnk_vals.prod().prod()
        testStats.append(testStat)
        count = rnk_vals.count().sum()
        counts.append(count)


    csR = ocl.dfCS.ix[brd][ind]



    cpRank = self.dfRank.ix[brd]
    cpSmRank = cpRank/100 # convert ranks back to 0 to 1
    nCPs.append(cpRes.shape[0])
    meanSer = cpRes.mean()
    meanRnk = cpRank.mean()
    prodRnk = cpSmRank.product()

for q in self.row_groups.keys():
    print "Connecting query group %s" % q
    for g in self.col_groups.keys():
        rnk_vals = self.query.pctrank.ix[self.row_groups[q],self.col_groups[g]]
        score_vals = self.query.score.ix[self.row_groups[q],self.col_groups[g]]





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



### oracle example 

file = '/xchip/cogs/rogerhu/scratch/test_oracle/result_SPEARMAN_n14x10.gctx'
# rp = connectors.RankProductConnector()
qres = queryresult.QueryResult()
qres.read(file)
ocl = oracle.Oracle(qres,out='/xchip/cogs/rogerhu/scratch/oracle')
ocl.build_groups(["pert_iname"], ["pert_iname"])
ocl.connect()
ocl.plot_all(connection_plots.rank_plot,variant="cell_id",plot_title="rank by cell_id")
ocl.mk_summary_html()
