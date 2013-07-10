 #! /usr/bin/env python

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd

work_dir = '/xchip/cogs/projects/target_id/CTD2_25June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

### make dictionary of BRDs (keys) and target gene symbols (values) - target_dict
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

### Use dgo to test drug-gene connectivity
reload(dgo)
targetDict = {'BRD-A35588707' : targetDict['BRD-A35588707'],
              'BRD-K06854232' : targetDict['BRD-K06854232']}
pDescDict = {'BRD-A35588707' : pDescDict['BRD-A35588707'],
             'BRD-K06854232' : pDescDict['BRD-K06854232']}

dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='KD',is_gold=True)
dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')
dg.test_known_connections(gp_type='KD',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='two_sided_dg_graphs_n10M',
                        make_graphs=True,
                        n_rand=1000,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='KD',
                metric='spearman',
                outName='apriori_two_sided_pass_FDR_n10M',
                alpha=0.2,
                make_graphs=True)
dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR_n10M')
dg.store_parameters_rpt()
outF = os.path.join(dg.outputdir,'drug-target_summary.txt')
dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')
