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


reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_connection')
dg.add_dictionary(targetDict=targetDictCGS)
# dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')

#make empty pDescDict
fullBRDs = []
for ind in dg.dfCS.index:
    brd = ind[0]
    fullBRDs.append(brd)
uniqBRDs = list(set(fullBRDs))
pDescDict = {}
for brd in uniqBRDs:
    pDescDict[brd] = '-666'

inFile = '/xchip/cogs/projects/oncoDome/OncoDome_genes.txt'
outDir = 'oncoDome_2July'
if not os.path.exists(work_dir+'/'+outDir):
    os.mkdir(work_dir+'/'+outDir)
cgsList = []
with open(inFile,'rt') as f:
    for string in f:
        splt = string[:-1]
        cgsList.append(splt)
geneAll = set(cgsList)
# check to see gene has a CGS
geneCGS = geneAll.copy()
for gene in geneAll:
	CM = mu.CMapMongo()
	CGSsigs = CM.find({'pert_type':'trt_sh.cgs','pert_iname':gene},{'sig_id':True,'pert_iname':True})
	if not CGSsigs:
		geneCGS.remove(gene)
for gene in geneCGS:
    dg.gene_to_drug_similarity(testGene=gene,
                                gp_type='KD',
                                metric='spearman',
                                outName=outDir + '/gene_to_drug',
                                pDescDict=pDescDict,
                                n_rand=10000,
                                n_uncorrected=20,
                                connection_test='two_sided')

