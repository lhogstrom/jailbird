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

### make target_dict
# targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_short.txt'
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

### test KD
reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')
# dg.test_known_connections(gp_type='KD',metric='spearman',pDescDict=pDescDict,make_graphs=True)
# dg.FDR_correction(pDescDict=pDescDict,metric='spearman',outName='apriori_connections_pass_FDR',alpha=0.2,make_graphs=False)
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
# dg.store_parameters_rpt()
# # dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # # dg.test_unknown_rank_product(gp_type='KD')
# # # # dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')
outF = os.path.join(dg.outputdir,'drug-target_summary_peyton.txt')
dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')

# #random
# geneList = ['ERBB2','MUC1','PIK3CA','MTOR','PPARG']
# geneList = ['MMP1']
# for gene in geneList:
# 	dg.gene_to_drug_similarity(testGene=gene,gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# #AURKA
# geneList = ['AURKA','AURKB','AURKAIP1']
inFile = '/xchip/cogs/projects/target_id/CTD2_25June2013/genes_with_connections.txt'
cgsList = []
with open(inFile,'rt') as f:
    for string in f:
        splt = string[:-1]
        cgsList.append(splt)
geneList = set(cgsList)
for gene in geneList:
    dg.gene_to_drug_similarity(testGene=gene,
                                gp_type='KD',
                                metric='spearman',
                                outName='gene_to_drug_connections',
                                pDescDict=pDescDict,
                                n_rand=10000,
                                n_uncorrected=20,
                                connection_test='two_sided')
# #drugbank connections
# geneList = ['PPARG,', 'FKBP1A', 'KIF11', 'MTOR', 'HMGCR', 'RRM1', 'ESR1', 'NR3C1', 'HMGCR', 'NNR3C1', 'HMGCR', 'NR3C1', 'PSMB1', 'PSMB5', 'RAF1', 'BRAF', 'CDK4', 'ESR1', 'NR3C1', 'NR3C1', 'NR3C1', 'NR3C1', 'R3C1', 'EGFR', 'HMGCR','EGFR', 'RRM1']


# # ### TEST OE
# # reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_OE_connection')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='OE',is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='OE',metric='spearman')
dg.test_known_connections(gp_type='OE',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='two_sided_dg_graphs',
                        n_rand=10000000,
                        make_graphs=False,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='OE',
                metric='spearman',
                outName='apriori_two_sided_pass_FDR_n10M',
                alpha=0.2,
                make_graphs=True,
                specificity_graph=True)
dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR_n10M',specificity_graph=True)
# dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # dg.test_unknown_rank_product(gp_type='KD')
# # # # dg.FDR_correction(pDe
outF = os.path.join(dg.outputdir,'drug-target_two_sided_summary.txt')
dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')
dg.store_parameters_rpt()


# # #re-asign varabile so computation is not lost when reloading the class
# # #dgCopy = dg
reload(dgo)
# dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
dg.dfRank = dgCopy.dfRank
dg.dfCS = dgCopy.dfCS
dg.pVec = dgCopy.pVec
dg.pDict = dgCopy.pDict
dg.countDict = dgCopy.countDict
dg.cellCountDict = dgCopy.cellCountDict
dg.nConnectionDict = dgCopy.nConnectionDict
dg.connectionsPassFDR = dgCopy.connectionsPassFDR
dg.pThreshFDR = dgCopy.pThreshFDR
dg.connection_test = dgCopy.connection_test

# # #hyperlink cmd
# outpath = '/'.join([dg.outputdir,'AURKA_gene_to_drug_connections'])
# hyperLnkPath = '/xchip/cogs/web/icmap/hogstrom/AURKA_KD_to_CTD2_connections'
# cmd = ' '.join(['ln -s',
# 		 outpath,
# 		 hyperLnkPath])
