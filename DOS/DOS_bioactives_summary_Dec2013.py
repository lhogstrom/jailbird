'''
-identify the bioactive DOS compounds in CMAP
-calculate summary stats
-prep for further analysis

Larson Hogstrom, 12/2013
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import cmap
import os

wkdir = '/xchip/cogs/projects/DOS/bioactivity_summary_Dec2013'

##grab anything that appeared on a DOS plate
# CM = mu.CMapMongo()
# dosQuery = CM.find({'sig_id':{'$regex':'DOS'},'pert_type':'trt_cp'}, #, 
#         {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
#         toDataFrame=True)
# potentialDos = set(dosQuery['pert_id'])

# #check the pert_collection status of each compounds
# mc = mu.MongoContainer()
# pertInfo = mc.pert_info.find({'pert_id':{'$in':list(potentialDos)}},
#             {},toDataFrame=True)
# collectionSets = set(pertInfo['pert_icollection'])


#get all cps from DOS collection
mc = mu.MongoContainer()
pertInfo = mc.pert_info.find({'pert_icollection':'DOS'},
            {},toDataFrame=True)
#check that it doesn't have a known pert_iname
inameSer = pertInfo['pert_iname']
inameFrm = pd.DataFrame(inameSer)
#which values do not start with BRD?
notBRDiname = pertInfo[~inameSer.str.contains('BRD')]
isBRDiname = pertInfo[inameSer.str.contains('BRD')]
dosBrds = isBRDiname['pert_id']

########################
## all DOS signatures ##
########################

# get signatures for all DOS compounds
CM = mu.CMapMongo()
dosQuery = CM.find({'pert_id':{'$in':list(dosBrds)},'pert_type':'trt_cp'}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
        toDataFrame=True)
dosSetLen = len(set(dosQuery['pert_id']))
dosGrped = dosQuery.groupby(['pert_id'])
countDict = {}
for grp in dosGrped:
    grpName = grp[0]
    cellSet = set(grp[1]['cell_id'])
    nCells = len(cellSet)
    countDict[grpName] = nCells
countSer = pd.Series(countDict)
countMax = max(countSer)
bins = np.arange(countMax+1)
plt.hist(countSer,bins)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('number of cell lines',fontweight='bold')
plt.xticks(bins)
plt.title('All DOS compounds (' + str(dosSetLen) + ') - cell lines tested tested')
outF = os.path.join(wkdir, 'DOS_cps_cell_line_distribution.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

#####################
## is_gold DOS ######
#####################

CM = mu.CMapMongo()
dosGold = CM.find({'pert_id':{'$in':list(dosBrds)},'pert_type':'trt_cp','is_gold':True}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
        toDataFrame=True)
dosGoldLen = len(set(dosGold['pert_id']))
dosGrped = dosGold.groupby(['pert_id'])
countGoldDict = {}
for grp in dosGrped:
    grpName = grp[0]
    cellSet = set(grp[1]['cell_id'])
    nCells = len(cellSet)
    countGoldDict[grpName] = nCells
countSerGold = pd.Series(countGoldDict)
countSerGold.name = 'n_cell_lines'
# goldCounts = dosGold.groupby(['pert_id','cell_id']).size().reset_index()
# countSerGold2 = goldCounts.groupby('pert_id').size()
countMax = max(countSerGold)
bins = np.arange(countMax+1)
plt.hist(countSerGold,bins)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('number of cell lines',fontweight='bold')
plt.title('gold DOS compounds (' + str(dosGoldLen) + ') - cell lines is_gold')
plt.xticks(bins)
outF = os.path.join(wkdir, 'DOS_isGold_cell_line_distribution.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

#which compounds are is_gold in more than one cell line?
multiGold = countSerGold[countSerGold > 1]

## which DOS compounds are in summly space?
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_sym_n7322x7322.gctx'
gt = gct.GCT()
gt.read(mtrxSummly)
summFrm = gt.frame
# find dos cps in summly space
summBrds = summFrm.index.values
summSet = set(summBrds)
dosSet = set(countSer.index)
overlapSet = dosSet.intersection(summSet)
#plot hist of cell counts - dos in summly 
overlapSer = countSerGold.reindex(list(overlapSet))
overlapCount = len(overlapSer)
countMax = max(overlapSer)
bins = np.arange(countMax+1)
plt.hist(overlapSer,bins)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('number of cell lines',fontweight='bold')
plt.title('summlySpace DOS compounds (' + str(overlapCount) + ') - cell lines is_gold')
plt.xticks(bins)
outF = os.path.join(wkdir, 'DOS_summly_cell_line_distribution.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

