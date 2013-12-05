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
columnPerts = gt.get_column_meta('pert_id')
summFrm = gt.frame
summFrm.columns = columnPerts
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

# ### for the 200 in summly space ---> run sig_introspect 

# ### run summly for bioactives - first in matched mode for the 234, then in independent mode

## how do gold signatures correlate across cell lines?
multiRes = dosGold[dosGold.pert_id.isin(multiGold.index)]
# load Expression Signatures
sigs = multiRes['sig_id'].values
afPath = cmap.score_path
gt = gct.GCT()
gt.read(src=afPath,cid=list(sigs),rid='lm_epsilon')
expFrm = gt.frame

pairCorr = np.corrcoef(expFrm.values.T)
corrFrm = pd.DataFrame(data=pairCorr,index=expFrm.columns,columns=expFrm.columns)

multiRes.index = multiRes['sig_id']
pertGrped = multiRes.groupby('pert_id')
#coherence of signatures across all cell lines
corrDict = {}
randDict = {}
observedCorrs = []
randCorrs = []
for brd in pertGrped.groups:
    sigs = pertGrped.groups[brd]
    smFrm = corrFrm.reindex(index=sigs,columns=sigs)
    nm = smFrm.shape[0]
    iUp = np.tril_indices(nm)
    smFrm.values[iUp] = np.nan
    corrArray = smFrm.unstack().values
    corrArray = corrArray[~np.isnan(corrArray)]
    corrDict[brd] = corrArray
    observedCorrs.extend(corrArray)
    # null distribution --> random pairs of sig_ids
    iRand = np.random.choice(multiRes.shape[0],nm,replace=False)
    randSigs = multiRes.index[iRand]
    rndFrm = corrFrm.reindex(index=randSigs,columns=randSigs)
    rndFrm.values[iUp] = np.nan
    rndArray = rndFrm.unstack().values
    rndArray = rndArray[~np.isnan(rndArray)]
    randDict[brd] = rndArray
    randCorrs.extend(rndArray)
# make histogram of self-self connections vs. random 
h1 = plt.hist(observedCorrs,30,color='b',range=[-1,1],label=['observed'],alpha=.4)
h2 = plt.hist(randCorrs,30,color='r',range=[-1,1],label='random',alpha=.4)
plt.legend()
plt.xlabel('corr')
plt.ylabel('freq')
plt.title('DOS compound self connections across cell lines')
outF = os.path.join(wkdir, 'DOS_pairwise_self_connection_distribution.png')
plt.savefig(outF, bbox_inches=0)
plt.close()

FileWtcs = '/xchip/cogs/projects/connectivity/query/wtcs.lm50/sim_wtcs.lm50_COMBINED.gctx'
gt = gct.GCT()
gt.read(src=FileWtcs,cid=list(sigs),rid=list(sigs))
wtcsFrm = gt.frame

# cluster the summly space for 
# dosSumm = summFrm.reindex(index=list(overlapSet),columns=list(overlapSet))

# plt.imshow(dosSumm.values,
#         interpolation='nearest',
#         cmap=matplotlib.cm.RdBu_r) #,
#         vmin=0, 
#         vmax=1)

# import scipy.cluster.hierarchy as hcluster

# seq_array = np.transpose(np.array(dosSumm.values))
# pdist = scipy.spatial.distance.pdist(seq_array)
# z = scipy.cluster.hierarchy.single(pdist)

# #compute the histogram of the linkage distances and find the first empty bin. Use that
# #bin to set the cutoff for cluster separation
# zdist = [z[i,2] for i in range(len(z))]
# n,bins,patches = plt.hist(zdist,bins=bins) #@UnusedVariable
# plt.clf()
# try:
#     distance_cutoff = bins[np.where(n == 0)][0]
# except IndexError:
#     distance_cutoff = np.max(bins)

# #cut the linkage hierarchy at the computed distance cutoff and assign cluster membership
# t = scipy.cluster.hierarchy.fcluster(z, distance_cutoff, criterion='distance')
# #build a 2D list in which each entry is a list of data points that fall in each cluster
# num_clusters = np.max(t)
# clusters = []
# cluster_data = []
# for i in range(num_clusters):
#     clusters.append([labels[j] for j in range(len(t)) if t[j] == i+1])
#     cluster_data.append([input_sequence[j] for j in range(len(t)) if t[j] == i+1])

# #return the list of clusters
# return (clusters,cluster_data)


