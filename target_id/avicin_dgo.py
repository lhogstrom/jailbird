import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd

work_dir = '/xchip/cogs/projects/avicins/dgo_result'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

avicinsBrds = ['BRD-A15100685','BRD-A33746814','BRD-A69592287','BRD-A70150975'] #avicin-d, avicin-g, oxetane, hydroxyl,
pDescDict = {'BRD-A15100685':'avicin-d','BRD-A33746814':'avicin-g','BRD-A69592287':'oxetane','BRD-A70150975':'hydroxyl'}

goiTested = ['PIK3CA',
     'PIK3CB',
     'GAPDH',
     'AKT1',
     'AKT2',
     'MTOR',
     'ALDOA',
     'NFKB1',
     'MYC']

# create target dictionary for genes of interest
targetDict = {}
for brd in avicinsBrds:
    targetDict[brd] = goiTested

### test KD
dg = dgo.QueryTargetAnalysis(out=work_dir + '/KD_spearman_all_doses')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='KD',is_gold=True)
dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')
dg.test_known_connections(gp_type='KD',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='test_dg_graphs2',
                        conn_thresh=.05,
                        make_graphs=True,
                        n_rand=100000,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='KD',
                metric='spearman',
                outName='test_FDR2',
                alpha=0.2,
                make_graphs=True,
                specificity_graph=True)
dg.fdr_html_summary(fdrDir='test_FDR2',specificity_graph=True)

## TEST OE
reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/OE_spearman_all_doses')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='OE',is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='OE',metric='spearman')
dg.test_known_connections(gp_type='OE',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='two_sided_dg_graphs',
                        n_rand=100000,
                        make_graphs=True,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='OE',
                metric='spearman',
                outName='apriori_two_sided_pass_FDR',
                alpha=0.2,
                make_graphs=True,
                specificity_graph=False)
dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR',specificity_graph=True)
