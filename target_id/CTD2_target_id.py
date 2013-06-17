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

work_dir = '/xchip/cogs/projects/target_id/CTD2_17June2013b'
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

test1 = 'OEB001_A375_96H:BRDN0000399163:-666' #set random sig_id to initialize dgo object
test2 = 'OEB001_A375_96H:BRDN0000400484:-666'
### test KD
# reload(dgo)
# dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_KD_connection')
# dg.add_dictionary(targetDict=targetDict)
# # dg.get_sig_ids(genomic_pert='KD')
# # dg.run_drug_gene_query(max_processes=10)
# # #wait until queries finish
# dg.make_result_frames(gp_type='KD')
# dg.test_known_connections(gp_type='KD',pDescDict=pDescDict)
# dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_apriori_connections')
# dg.test_unknown_rank_product(gp_type='KD')
# dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')

### TEST OE
reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_OE_connection')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='OE')
dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='OE')
dg.test_known_OE_connections(pDescDict=pDescDict,gp_type='OE')
dg.FDR_correction(pDescDict=pDescDict)


# ### debug creation of OE connection graphs

#     def test_known_OE_connections(self,pDescDict,gp_type='OE',n_rand=100000,make_graphs=True):
#         '''
#         test known connections saved in targetDict across cell lines
#         --> retrieve similarity scores
#         --> calculate p-values (based on percent rank)
#         --> generate graphs of connections

#         inputs:
#         pDescDict = dictionary of common names for brds tested
#         n_rand = number of permutations to calculate null distribution of test statistic
        
#         outputs:
#         pDict = dictionary of p-values for each drug-gene relationship tested
#         pVec = list of p-values for each drug-gene relationship tested
#         countDict = dictionary of instances per drug-gene relationship
#         brdSkipped = compounds skipped due to lack of data
#         cgsSkipped = genomic perturbations skipped due to lack of data
#         '''
#         #set output dir
#         gp_type='OE'
#         n_rand=10000
#         make_graphs=True
#         graphDir = dg.outputdir + '/drug_target_graphs'
#         if not os.path.exists(graphDir):
#             os.mkdir(graphDir)
#         #get brds from result dataframe
#         brdSkipped = []
#         cgsSkipped = []
#         BRDsTested = []
#         for ind in dg.dfRank.index:
#             brd = ind[0]
#             BRDsTested.append(brd)
#         brdRsltSet = set(BRDsTested)
#         #get cgs tested 
#         # cols = dg.dfRank.columns
#         cols = []
#         for col in dg.dfRank.columns:
#             cols.append(col[0])        
#         countDict = {}
#         pDict = {}
#         pVec = []
#         prog = progress.DeterminateProgressBar('Connection test')
#         for ibrd,brd in enumerate(dg.targetDict):
#             # skip pert if not in result file
#             prog.update(brd,ibrd,len(dg.targetDict))
#             if not brd in brdRsltSet:
#                 brdSkipped.append(brd)
#                 continue
#             targets = dg.targetDict[brd]
#             cpRes = dg.dfCS.ix[brd]
#             cpRank = dg.dfRank.ix[brd]
#             meanSer = cpRes.mean()
#             meanRnk = cpRank.mean()
#             nullCnt = pd.isnull(cpRes)
#             #how many cell lines were both the pert and target tested in
#             valCounts = nullCnt.shape[0] - nullCnt.sum(axis=0)
#             for target in targets:
#                 tarList = [inst for inst in cols if inst.split('_')[0] == target]
#                 if len(tarList) == 0: #skip if drug target not tested
#                     cgsSkipped.append(target)
#                     continue
#                 for ind in tarList:
#                     rnkSer = cpRank[ind]
#                     rnkSer = rnkSer.unstack()
#                     rnkSer = rnkSer[rnkSer.notnull()]
#                     csSer = cpRes[ind]
#                     csSer = csSer.unstack()
#                     csSer = csSer[csSer.notnull()]
#                     #skip if cgs not tested in the same cell line as cp
#                     if len(rnkSer) == 0:
#                         cgsSkipped.append(ind)
#                         continue
#                     ### calculate p-value - based on percent rank products
#                     rnkSmll = rnkSer/100
#                     testStat = rnkSmll.prod()
#                     n_obs = rnkSer.shape[0]
#                     # theoretical null
#                     ### simulate random draws from percent rank list
#                     permMtrx = np.random.rand(n_obs,n_rand)
#                     nullDist = permMtrx.prod(axis=0)
#                     #number of null values more extreme than observed (one sided)
#                     exVals = nullDist[nullDist<testStat]
#                     nExtreme = len(exVals)
#                     pVal = (nExtreme+1)/float(len(nullDist))
#                     pVec.append(pVal)
#                     pDict[brd + '-' + ind] = pVal
#                     #make summary output
#                     outF = os.path.join(graphDir,brd +'_' + ind + '_drug-target_summary.txt')
#                     # dg.__make_CS_summary(brd,pDescDict[brd],rnkSer,csSer,pVal,outF,gp_type)
#                     if make_graphs:
#                         ### cs wadden gram
#                         sKeysStr = []
#                         count = 0
#                         for i,cs in enumerate(csSer):
#                             if pd.isnull(cs):
#                                 continue
#                             else:
#                                 count = count + 1
#                                 sKeysStr.append(csSer.index[i][1].split('_')[1])
#                                 yVals = count
#                                 plt.scatter(cs,yVals)
#                         plt.xlim((-1, 1))
#                         plt.ylim((0,count+1))
#                         plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
#                         plt.xlabel('wtcs')
#                         plt.ylabel('cell line')
#                         plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type)
#                         plt.savefig(os.path.join(graphDir,brd +'_' + ind + '_connections.png'))
#                         plt.close()
#                         #rank wadden gram
#                         sKeysStr = []
#                         count = 0
#                         rnkList = cpRank[ind]
#                         for i,rnk in enumerate(rnkSer):
#                             if pd.isnull(rnk):
#                                 continue
#                             else:
#                                 count = count + 1
#                                 sKeysStr.append(csSer.index[i][1].split('_')[1])
#                                 yVals = count
#                                 plt.scatter(rnk,yVals)
#                         plt.xlim((0, 100))
#                         plt.ylim((0,count+1))
#                         plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
#                         plt.xlabel('percent rank')
#                         plt.ylabel('cell line')
#                         plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type)
#                         plt.savefig(os.path.join(graphDir,brd +'_' + ind + '_percent_rank.png'))
#                         plt.close()
#         dg.brdSkipped = brdSkipped
#         dg.cgsSkipped = cgsSkipped
#         dg.pDict = pDict
#         dg.pVec = pVec


