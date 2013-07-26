#! /usr/bin/env python

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo_oracle as dgo
import cmap.util.progress as progress
import pandas as pd

work_dir = '/xchip/cogs/projects/target_id/oracle_16Jul'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

pert_list = ['BRD-K69840642', 'BRD-K41859756']
pDescDict = {'BRD-K69840642':'ISOX','BRD-K41859756':'NVP-AUY922'}
targetDict = {'BRD-K69840642':['HSPA9'],'BRD-K41859756':['HSPA5']}

# dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
# dg = dgo.make_qres_obj(metric='wtcs',gp_type='KD')
# dg.get_sig_ids(genomic_pert='KD',targetDict_loaded=False,pert_list=pert_list,is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)

### test KD
reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')
dg.make_qres_obj(metric='spearman',gp_type='KD')
# dg.test_known_connections(gp_type='KD',metric='spearman',pDescDict=pDescDict,make_graphs=True)
# dg.FDR_correction(pDescDict=pDescDict,metric='spearman',outName='apriori_connections_pass_FDR',alpha=0.2,make_graphs=False)

# dg.test_known_connections_old(gp_type='KD',
#                         metric='spearman',
#                         pDescDict=pDescDict,
#                         outName='apriori_connection_old',
#                         conn_thresh=.05,
#                         make_graphs=True,
#                         n_rand=10000,
#                         connection_test='two_sided')
dg.test_known_connections(gp_type='KD',
                        metric='spearman',
                        outName='apriori_connection_graphs',
                        conn_thresh=.05,
                        make_graphs=True,
                        niter=10000,
                        connection_direction='both')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='KD',
                metric='spearman',
                outName='test_FDR2',
                alpha=0.2,
                make_graphs=True,
                specificity_graph=True)
dg.fdr_html_summary(fdrDir='test_FDR2',specificity_graph=True)

# dgCopy = dg
reload(dgo)
dg.qres = dgCopy.qres


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