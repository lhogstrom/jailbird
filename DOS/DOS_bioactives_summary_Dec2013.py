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
from matplotlib.ticker import NullFormatter

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
dosQuery.index = dosQuery['sig_id']
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
dosGold.index = dosGold['sig_id']
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

########################
## make SC Plots #######
########################

x = dosGold['distil_cc_q75'].values
y = dosGold['distil_ss'].values
# y = nsample
nullfmt   = NullFormatter()         # no labels
# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left+width+0.02
#
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]
# start with a rectangular Figure
plt.figure(1, figsize=(8,8))
axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)
# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)
# the scatter plot:
axScatter.scatter(x, y, marker=".", alpha=.1) #,marker='.'
binwidth = 0.25
# xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
xymax = np.max(np.fabs(y))
# lim = ( int(xymax/binwidth) + 1) * binwidth
lim = ( int(xymax/binwidth) + xymax) * binwidth
#set 
xMin = -.3
xMax = 1
xBinwidth = .025
axScatter.set_xlim( (xMin, xMax) )
axScatter.set_ylim( (0, lim) )
#make hist
binsX = np.arange(xMin, xMax + xBinwidth, xBinwidth)
binsY = np.arange(0, lim + binwidth, binwidth)
# axHistx.hist(x, bins=bins)
axHistx.hist(x, bins=binsX, range=[xMin,xMax])
axHisty.hist(y, bins=binsY, orientation='horizontal')
axHistx.axis('off')
axHisty.axis('off')
#set lims
axScatter.set_ylabel('ss',fontweight='bold')
axScatter.set_xlabel('CC_q75',fontweight='bold')
# axScatter.set_title('DMSO signatures - LINCS core cell lines (n = 6471)')
outF = os.path.join(wkdir, 'DOS_SC_plot_all_is_gold.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()




## which DOS compounds are in summly space?
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_sym_n7322x7322.gctx'
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_n7147x7147.gctx'
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_lass_n39560x7147.gctx'
gt = gct.GCT()
# gt.read_gctx_col_meta(mtrxSummly)
# gt.read_gctx_row_meta(mtrxSummly)
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

#which compounds are in independent mode space?
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_lass_n39560x7147.gctx'
gt = gct.GCT()
gt.read_gctx_col_meta(mtrxSummly)
gt.read_gctx_row_meta(mtrxSummly)
indSummSigs = gt.get_column_meta('sig_id')
indSummInames = gt.get_column_meta('pert_iname')
sigSer = pd.Series(index=indSummSigs, data=indSummInames)
#which dos cps are in the summly indpend
summBrds = set(sigSer.index)
goldDosBrds = set(dosGold['sig_id'].values)
summGold = summBrds.intersection(goldDosBrds)


dosSer = sigSer.reindex(dosGold['sig_id'].values)


### make summary table:
# 1) pert_id
# 2) times_profiled_in_a2
# 3) times_gold_in_a2
# 4) is_gold_cell lines
goldGrped = dosGold.groupby('pert_id')
allGrped = dosQuery.groupby('pert_id')
#with vcap
sFrm = pd.DataFrame()
for brd in multiGold.index:
    allSigs = allGrped.groups[brd]
    nTested = len(allSigs)
    goldSigs = goldGrped.groups[brd]
    nGold = len(goldSigs)
    goldCells = dosGold.ix[goldSigs,'cell_id']
    goldCellSet = set(goldCells)
    noVCAP = goldCellSet.copy()
    if 'VCAP' in noVCAP:
        noVCAP.remove('VCAP')
    #which compounds are is gold in two cell lines OTHER than VCAP
    if len(noVCAP) > 1:
        nDict = {}
        nDict['times_profiled_in_a2'] = nTested
        nDict['times_gold_in_a2'] = nGold
        nDict['is_gold_cell_lines'] = list(goldCellSet)
        brdFrm = pd.DataFrame({brd:nDict})
        sFrm = pd.concat([sFrm,brdFrm],axis=1)
sFrm = sFrm.T
#convert list of cell lines to a comma seperated strings
sFrm['is_gold_cell_lines'] = sFrm['is_gold_cell_lines'].str.join(',')
sFrm = sFrm.reindex(columns=['times_profiled_in_a2','times_gold_in_a2','is_gold_cell_lines'])
outF = os.path.join(wkdir, 'DOS_gold_cell_line_summary.txt')
sFrm.to_csv(outF,sep='\t', header=True)

### for the 200 in summly space ---> run sig_introspect 
outIntrospect = wkdir + '/sig_introspect_DOS_summly_space_non_gold'
if not os.path.exists(outIntrospect):
    os.mkdir(outIntrospect)
#get all sig_ids for dos cps in summlyspace
CM = mu.CMapMongo()
summSpaceQuery = CM.find({'pert_id':{'$in':list(overlapSet)},'pert_type':'trt_cp'}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
        toDataFrame=True)
qSer = summSpaceQuery['sig_id']
outF = outIntrospect + '/sig_ids_dos_summSpace.grp'
qSer.to_csv(outF,index=False,header=False)
#run sig_introspect
# cmd = ' '.join(['rum -q hour',
#      '-d sulfur_io=100',
cmd = ' '.join(['rum -q hour',
     '-o ' + outIntrospect,
     '-x sig_introspect_tool ',
     '--sig_id ' + outF,
     '--query_group pert_id',
     '--metric wtcs',
     '--out ' + outIntrospect])
os.system(cmd)
#which compounds have a pre-calculated introspect result?
resPreCalc = '/xchip/cogs/projects/connectivity/introspect/'

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


