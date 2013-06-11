#! /usr/bin/env python

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo

work_dir = '/xchip/cogs/projects/target_id/5June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)
### make target_dict
# targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_short.txt'
targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_drug_targets.txt'
'/xchip/cogs/projects/target_id/7June2014/A2_DrugBank_targets_tab.txt'
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
			if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
				continue
			else:
				targetDict[pID] = targets
				pDescDict[pID] = pDesc

db = mu.CMapMongo()
test1 = db.find({'cell_id':'A375','is_gold':True,'pert_type':'trt_oe'},{'sig_id':1})
test2 = db.find({'cell_id':'A375','is_gold':True,'pert_type':'trt_sh'},{'sig_id':1})
test1 = test1[1:10]
test2 = test2[1:10]
# t = dgo.Oracle(test1,test2,out=work_dir + '/Oracle')
# t.compute_scores()
# t.get_results()

reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_gene')
dg.add_dictionary(targetDict=targetDict)
# dg.get_drug_kd_sig_ids()
# dg.run_drug_gene_query()
dg.make_result_frames()
dg.test_known_connections(pDescDict=pDescDict,make_graphs=False)
dg.FDR_correction(pDescDict=pDescDict)




### test rand vec time
# n_obs = 44
# n_rand = 1000000
# permMtrx = np.random.rand(n_obs,n_rand)
# nullDist = permMtrx.prod(axis=0)
# testStat = nullDist[4]
# exVals = nullDist[nullDist<testStat]
# plt.hist(nullDist,bins=np.logspace(-18, 0, 50))
# plt.gca().set_xscale('log')
# plt.show()

