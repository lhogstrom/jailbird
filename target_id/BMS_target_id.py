#! /usr/bin/env python
'''
analyze drug connections with BMS cps
'''

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress

work_dir = '/xchip/cogs/projects/target_id/BMS_11June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

# Query instances of BMS cps
CM = mu.CMapMongo()
pert_Q = CM.find({'pert_iname':{'$regex':'BMS'},'pert_type':'trt_cp'},{'sig_id':True,'pert_iname':True,'pert_id':True,})

#make brd-iname dictionary
pDescDict = {}
for sig in pert_Q:
	pDescDict[sig['pert_id']] = sig['pert_iname']
pert_list = pDescDict.keys()

### run dgo object
test1 = 'OEB001_A375_96H:BRDN0000399163:-666' #set random sig_id to initialize dgo object
test2 = 'OEB001_A375_96H:BRDN0000400484:-666'
reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_KD_connection')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD',targetDict_loaded=False,pert_list=pert_list)
# dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
dg.make_result_frames()
dg.test_unknown_connections(gp_type='KD',pDescDict=pDescDict)
# dg.FDR_correction(pDescDict=pDescDict)

### run OE analysis
reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_OE_connection')
# dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='OE',targetDict_loaded=False,pert_list=pert_list)
dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='OE')
dg.test_unknown_connections(gp_type='OE',pDescDict=pDescDict)
# dg.FDR_correction(pDescDict=pDescDict)
