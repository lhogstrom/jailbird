 #! /usr/bin/env python

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo

work_dir = '/xchip/cogs/projects/PRISM/target_id/12June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

#load in gene/ cp list and make target dict
geneF = '/xchip/cogs/projects/PRISM/target_id/gene_targets.txt'
geneList = [line.strip() for line in open(geneF)]
## make target dict 
pDescDict = { "BRD-K59369769": "VX-680", "BRD-K36740062": "GSK-1070916","CMAP-AZD-1152HQPA": "AZD-1152HQPA","BRD-K29830875": "-666","BRD-K01737880": "-666","BRD-K83963101": "MLN8054"}
targetDict = {}
for pert in pDescDict:
	targetDict[pert] = geneList

test1 = 'OEB001_A375_96H:BRDN0000399163:-666' #set random sig_id to initialize dgo object
test2 = 'OEB001_A375_96H:BRDN0000400484:-666'
### test KD
reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_KD_connection')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD')
# dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD')
# dg.test_known_connections(gp_type='KD',pDescDict=pDescDict)
dg.test_unknown_rank_product(gp_type='KD')
dg.FDR_correction(pDescDict=pDescDict)

### TEST OE
# dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_OE_connection')
# dg.add_dictionary(targetDict=targetDictCGS)
# dg.get_sig_ids(genomic_pert='OE')
# # dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
# dg.make_result_frames(gp_type='OE')
# dg.test_known_OE_connections(pDescDict=pDescDict,gp_type='OE')
# dg.FDR_correction(pDescDict=pDescDict)