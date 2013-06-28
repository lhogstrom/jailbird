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

work_dir = '/xchip/cogs/projects/target_id/24June_Barzily_Rokni_DOS_full'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

### make target_dict
targetSheetF = '/xchip/cogs/projects/target_id/Barzily_Rokni_DOS/BR_drug_target_list.txt'
# targetSheetF = '/xchip/cogs/projects/target_id/Barzily_Rokni_DOS/BR_drug_target_list_short.txt'
targetDict = {}
pDescDict = {}
with open(targetSheetF,'rt') as f:
	for string in f:
		splt = string.split('\r')
		for i,line in enumerate(splt):
			splt2 = line.split('\t')
			pID = splt2[0][:13] #the pert_id listed the line
			pDesc = splt2[1]
			targets = splt2[2]
			targets = targets.split(';')
			targets = [x[:-1] for x in targets if x != '']
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
dg.test_known_connections(gp_type='KD',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='two_sided_dg_graphs',
                        make_graphs=False,
                        n_rand=1000000,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='KD',
                metric='spearman',
                outName='apriori_two_sided_pass_FDR',
                alpha=0.2,
                make_graphs=True)
# dg.FDR_correction(pDescDict=pDescDict,metric='spearman',outName='apriori_connections_pass_FDR',alpha=0.2,make_graphs=True)
dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR')
outF = os.path.join(dg.outputdir,'drug-target_summary_two_sided.txt')
dg.make_target_summary(outF,dir_loc='apriori_two_sided_pass_FDR')

# dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # dg.test_unknown_rank_product(gp_type='KD')
# # # dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')

# #random
# geneList = ['ERBB2','MUC1','PIK3CA','MTOR','PPARG']
# for gene in geneList:
# 	dg.gene_to_drug_similarity(testGene=gene,gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# #AURKA
# geneList = ['AURKA','AURKB','AURKAIP1']
# for gene in geneList:
# 	dg.gene_to_drug_similarity(testGene=gene,gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# #drugbank connections
# geneList = ['PPARG,', 'FKBP1A', 'KIF11', 'MTOR', 'HMGCR', 'RRM1', 'ESR1', 'NR3C1', 'HMGCR', 'NNR3C1', 'HMGCR', 'NR3C1', 'PSMB1', 'PSMB5', 'RAF1', 'BRAF', 'CDK4', 'ESR1', 'NR3C1', 'NR3C1', 'NR3C1', 'NR3C1', 'R3C1', 'EGFR', 'HMGCR','EGFR', 'RRM1']


# # ### TEST OE
# # reload(dgo)
# # dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_OE_connection')
# # dg.add_dictionary(targetDict=targetDict)
# # dg.get_sig_ids(genomic_pert='OE')
# # # dg.run_drug_gene_query(max_processes=10)
# # # #wait until queries finish
# # dg.make_result_frames(gp_type='OE')
# # dg.test_known_OE_connections(pDescDict=pDescDict,gp_type='OE',make_graphs=False)
# # dg.FDR_correction(pDescDict=pDescDict,gp_type='OE',outName='FDR_pass_graphs_alpha2',alpha=0.2)
# # dg.fdr_html_summary(fdrDir='FDR_pass_graphs_alpha2')

# # #re-asign varabile so computation is not lost when reloading the class
# # #dgCopy = dg
reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
dg.dfRank = dgCopy.dfRank
dg.dfCS = dgCopy.dfCS
dg.pVec = dgCopy.pVec
dg.pDict = dgCopy.pDict
dg.connectionsPassFDR = dgCopy.connectionsPassFDR
dg.cellCountDict = dgCopy.cellCountDict
dg.countDict = dgCopy.countDict
dg.nConnectionDict = dgCopy.nConnectionDict
# # #hyperlink cmd
# outpath = '/'.join([dg.outputdir,'AURKA_gene_to_drug_connections'])
# hyperLnkPath = '/xchip/cogs/web/icmap/hogstrom/AURKA_KD_to_CTD2_connections'
# cmd = ' '.join(['ln -s',
# 		 outpath,
# 		 hyperLnkPath])


### make a list of cell lines in which the cps were tested
outfile = work_dir + '/drug_cell_lines_tested.txt'
headers = ['brd','cell_lines']
with open(outfile,'w') as f:
	f.write('\t'.join(headers) + '\n')
	prog = progress.DeterminateProgressBar('perturbation cid query')
	for ibrd,brd in enumerate(targetDict):
		prog.update('querying cps',ibrd,len(targetDict))
		CM = mu.CMapMongo()
		#is gold
		# pert_query = CM.find({'pert_id':{'$regex':brd},'is_gold':True},{'sig_id':True,'cell_id':True,'pert_id':True,'is_gold':True})
		# also return non is_golds
		pert_query = CM.find({'pert_id':{'$regex':brd}},{'sig_id':True,'cell_id':True,'pert_id':True,'is_gold':True})
		if pert_query:
			cells = [q['cell_id'] for q in pert_query]
			cells.insert(0,brd)
			f.write('\t'.join(cells) + '\n')
		else:
			f.write(brd + '\n')



