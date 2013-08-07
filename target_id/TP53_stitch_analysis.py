#! /usr/bin/env python
'''
Acording to STITCH, what are the genes that relate to the to the TP53 drugs?

how do the summly results relate to these?

'''
import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd
import matplotlib

work_dir = '/xchip/cogs/hogstrom/analysis/summly/TP53_dgo'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

stitchF = '/xchip/cogs/hogstrom/notes/stitch_db/STITCH_drug_target_800%2B_activators.txt'
sF = pd.read_csv(stitchF,delimiter='\t') # stitch frame

#generate list of drugs associated with p53
geneF = sF[sF['GeneSymbol'] == 'TP53']
pIDs = list(geneF['pert_id'].values)
pIDs = list(set(pIDs))
targetDict = {}
for brd in pIDs:
    targetDict[brd] = 'TP53'

#make dictionary of inames
pDescDict = {}
CM = mu.CMapMongo()
for pert in pIDs:
        pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True},{'pert_iname':True},limit=1)
        if pert_query:
            pDescDict[pert] = pert_query[0]
        else:
            pDescDict[pert] = '-666'


#run dgo - KDs
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='KD',is_gold=True)
dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')
# dg.test_known_connections(gp_type='KD',metric='spearman',pDescDict=pDescDict,make_graphs=True)
# dg.FDR_correction(pDescDict=pDescDict,metric='spearman',outName='apriori_connections_pass_FDR',alpha=0.2,make_graphs=False)
dg.test_known_connections(gp_type='KD',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='dg_graphs',
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
dg.store_parameters_rpt()
outF = os.path.join(dg.outputdir,'drug-target_summary_peyton.txt')
dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')

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
                        n_rand=100000,
                        make_graphs=False,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='OE',
                metric='spearman',
                outName='apriori_two_sided_pass_FDR',
                alpha=0.2,
                make_graphs=True,
                specificity_graph=True)
dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR',specificity_graph=True)
# dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # dg.test_unknown_rank_product(gp_type='KD')
# # # # dg.FDR_correction(pDe
outF = os.path.join(dg.outputdir,'drug-target_two_sided_summary.txt')
dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')
dg.store_parameters_rpt()


### run queries for cps
cgsCells = ['A375', 'A549', 'ASC', 'HA1E', 'HEPG2', 'HCC515', 'HT29', 'MCF7', 'NPC', 'PC3', 'VCAP']
processes = set()
file1 = work_dir + '/' + 'TP53_stitch_cp_sigs.grp'
for pert in pIDs:
    sigs = []
    for cell in cgsCells:
        CM = mu.CMapMongo()
        #dose min
        # pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell,'pert_dose':{'$gt':5}},{'sig_id':True},limit=1)
        #no dose min
        pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell},{'sig_id':True},limit=1)
        if pert_query:
            sigs.append(pert_query[0])
        # else:
        #     print 'no sig for ' + pert + ' in ' + cell
    with open(file1, 'a') as f:
        [f.write(x + '\n') for x in sigs]
metric = 'wtcs'
queryDir = work_dir + '/ctd2_sig_query'
if not os.path.exists(queryDir):
    os.mkdir(queryDir)
cmd = ' '.join(['rum -q local -f sig_query_tool',
         '--sig_id ' + file1,
         '--metric ' + metric,
         '--column_space full',
         '--out ' + queryDir,
         '--mkdir false',
         '--save_tail false'])
         # '--row_space bing', 
os.system(cmd)

