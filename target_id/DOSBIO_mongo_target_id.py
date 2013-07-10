
#! /usr/bin/env python
'''
analyze the DOSBIO plates - combine with data in cmap database to perform query

use DOS signatures generate queries of the CGS data (cell line specific results)
'''

import os
import cmap.io.gct as gct
import cmap.analytics.sc as sc
import glob as glob
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import cmap.analytics.dgo as dgo

work_dir = '/xchip/cogs/projects/DOS/8July_target_id'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

### get DOS BIO cps from mongo
CM = mutil.CMapMongo()
# pert_List = CM.find({'sig_id':{'$regex':'DOSBIO'},'pert_iname':{'$regex':'BRD'}},{'sig_id':True,'pert_id':True,'pert_iname':True})
pert_List = CM.find({'sig_id':{'$regex':'DOSBIO'},'pert_iname':{'$regex':'BRD'}},{'pert_id':True})
dosbioSet = set(pert_List)
# check to make sure the brds are DOS compounds and don't represent known compounds
inameDict = {}
for brd in dosbioSet:
    inames = CM.find({'pert_id':brd},{'pert_iname':True})
    inameSet = set(inames)
    inameDict[brd] = inameSet

### test KD
reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
# dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='KD',
            targetDict_loaded=False,
            pert_list=list(dosbioSet),
            is_gold=True)
dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')







