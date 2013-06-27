#! /usr/bin/env python
'''
analyze drug connections with erbb2 knockdown
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

work_dir = '/xchip/cogs/projects/target_id/ERBB2_5June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

# Query instances of cps
# CM = mu.CMapMongo()
# # pert_List = CM.find({'pert_iname':{'$regex':'5240'}},{'sig_id':True,'pert_iname':True,'pert_id':True,})
# # erbb2Lst = CM.find({'pert_iname':{'$regex':'ERBB2'},'pert_type':'trt_sh.cgs'},{'sig_id':True,'pert_iname':True,'pert_id':True,})
# erbb2Lst = CM.find({'pert_iname':{'$regex':'ERBB2'},'pert_type':'trt_oe'},{'sig_id':True,'pert_iname':True,'pert_id':True,})

#cps to check:
# AZD-8055 - BRD-K69932463
# AS-605240 - BRD-K41895714
# ASG05240
targetDict = {}
targetDict['BRD-K69932463'] = ['ERBB2']
targetDict['BRD-K41895714'] = ['ERBB2']

pDescDict = {}
pDescDict['BRD-K69932463'] = 'AZD-8055'
pDescDict['BRD-K41895714'] = 'AS-605240'

## test OE
reload(dgo)
dg = dgo.QueryTargetAnalysis(work_dir + '/drug_OE_connection')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='KD')
dg.run_drug_gene_query(max_processes=10)
#wait until queries finish
dg.make_result_frames(gp_type='KD')
dg.test_known_connections(pDescDict=pDescDict,gp_type='OE')
dg.FDR_correction(pDescDict=pDescDict)

### test KD
reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_KD_connection')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD')
# dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD')
# dg.test_known_connections(gp_type='KD',pDescDict=pDescDict)
# dg.FDR_correction(pDescDict=pDescDict)

## # do an umbiased search to see which drugs in CMAP best connect to HER2 CGS


