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
import os
import cmap.analytics.summly_null as SN
from statsmodels.distributions import ECDF
import cmap.io.gmt as gmt
import cmap.util.progress as update
from matplotlib import cm
import cmap.plot.colors as ccol
import scipy.cluster

wkdir = '/xchip/cogs/projects/DOS/PCL_comparison_Feb142014'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

########################
## get DOS collection ##
########################

def get_dos_BRDs():
    'return Series of all DOS compounds - use the DOS icollection'
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
    return dosBrds

########################
## all DOS signatures ##
########################

def get_dos_signatures(dosBrds):
    "1) return signature info for all DOS compounds \
    2) return list counts of number of cell lines tested"
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
    return dosQuery, countSer

############################
## summly matched DOS ######
############################

#which compounds are in independent mode space?
def get_summly_dos_indeces(dosBrds,mtrxSummly):
    "1) return dos compounds that are in summly matched space "
    gt = gct.GCT()
    gt.read(mtrxSummly)
    # indSummSigs = gt.get_column_meta('sig_id')
    # indSummInames = gt.get_column_meta('pert_iname')
    summFrm = gt.frame
    sigSer = pd.Series(index=summFrm.index, data=summFrm.columns)
    dosSer = sigSer[sigSer.index.isin(dosBrds)]
    return dosSer

##########################################
### run functions to retrieve DOS info ###
##########################################

dosBrds = get_dos_BRDs()
# dosQuery, countSer = get_dos_signatures(dosBrds)
# dosGold, countSerGold = get_dos_gold_signatures(dosBrds)
### use lass matrix
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_lass_n39560x7147.gctx'
# matrixType = 'rnkpt_indp_lass'
### use mrp4 mtrx
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_mrp4_n39560x7147.gctx'
# matrixType = 'mrp4'
### use lass matched matrix
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_n7147x7147.gctx'
matrixType = 'rnkpt_matched_lass'
dosSer = get_summly_dos_indeces(dosBrds,mtrxSummly)
outGRP = wkdir + '/dos_summly_ids.grp'
dosSer.to_csv(outGRP,index=False)
# nullMtrx = '/xchip/cogs/projects/connectivity/null/dmso/lass_n1000x7147.gctx'
nullMtrx = '/xchip/cogs/projects/connectivity/null/random/lass_n1000x7147.gctx'

##########################################
### get brds for the whole matrix ###
##########################################

gt = gct.GCT()
gt.read(mtrxSummly)
summFrm = gt.frame
pInames = gt.get_column_meta('pert_iname')
pIDs = gt.get_column_meta('pert_id')
pType = gt.get_column_meta('pert_type')
anntFrm = pd.DataFrame({'pert_id':pIDs,'pert_type':pType,'pert_iname':pInames},index=pIDs)
sigSer = pd.Series(index=summFrm.index, data=summFrm.columns)
outGRP = wkdir + '/summly_matched_ids.grp'
sigSer.to_csv(outGRP,index=False)

##########################################
### construct sig_cliquescore_tool command ###
##########################################

# groupGMT = '/xchip/cogs/projects/pharm_class/rnwork/cliques/cpd_groups_n147.gmt'
groupGMT = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140213/cliques.gmt'
cliqueGMT = gmt.read(groupGMT)
cliqFrm = pd.DataFrame(cliqueGMT)
cliqFrm['group_size'] = cliqFrm.sig.apply(len)
cliqFrm.index = cliqFrm['desc']

#observed 
cmd = ', '.join(['sig_cliquescore_tool(\'score\', \'' + mtrxSummly+ '\'',
        '\'summly_id\', \'' + outGRP + '\'',
         '\'clique\', \'' + groupGMT + '\'',
         '\'out\', \'' + wkdir + '\')'])

#null 
cmd = ', '.join(['sig_cliquescore_tool(\'score\', \'' + nullMtrx + '\'',
         '\'clique\', \'' + groupGMT + '\'',
         '\'out\', \'' + wkdir + '\')'])

#########################################
### load sig_cliquescore_tool results ###
#########################################

# cliques against DOS compounds
# cFile = '/xchip/cogs/projects/DOS/PCL_comparison_Feb162014/feb13/my_analysis.sig_cliquescore_tool.2014021313530491/clique_median_n145x234.gctx'
cFile = '/xchip/cogs/projects/DOS/PCL_comparison_Feb142014/feb14/my_analysis.sig_cliquescore_tool.2014021415040791/clique_median_n80x234.gctx'
gt = gct.GCT()
gt.read(cFile)
cliqDos = gt.frame
cliqDos = cliqDos.reindex(dosSer.values)
cliqDos.index = dosSer.index

# cliques against all compounds
# cFileFull = '/xchip/cogs/projects/DOS/PCL_comparison_Feb162014/feb14/my_analysis.sig_cliquescore_tool.2014021411101091/clique_median_n145x7147.gctx'
cFileFull = '/xchip/cogs/projects/DOS/PCL_comparison_Feb142014/feb14/my_analysis.sig_cliquescore_tool.2014021414592391/clique_median_n80x7147.gctx'
gt2 = gct.GCT()
gt2.read(cFileFull)
cliqFull = gt2.frame
cliqFull = cliqFull.reindex(sigSer.values)
cliqFull.index = sigSer.index
# transpose matrix and write
# cT = cliqFull.T 
# gt3 = gct.GCT()
# gt3.build_from_DataFrame(cT)

# cliques against DMSO
cFile = '/xchip/cogs/projects/DOS/PCL_comparison_Feb142014/feb14/my_analysis.sig_cliquescore_tool.2014021415010091/clique_median_n80x1000.gctx'
gt4 = gct.GCT()
gt4.read(cFile)
cliqDMSO = gt4.frame

# cliques against random signatures
cFile = '/xchip/cogs/projects/DOS/PCL_comparison_Feb142014/feb14/my_analysis.sig_cliquescore_tool.2014021415003091/clique_median_n80x1000.gctx'
gt5 = gct.GCT()
gt5.read(cFile)
cliqRnd = gt5.frame

#####################################
### plot DOS heatmap - make hist  ###
#####################################

cUnstack = cliqDos.unstack()
plt.hist(cUnstack,30)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('median ' + matrixType,fontweight='bold')
plt.title('distribution of DOS-clique connections')
outF = os.path.join(wkdir, 'dos_clique_connection_distribution.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

### heatmap 
ccol.set_color_map()
plt.imshow(cliqDos.values,
    interpolation='nearest',
    aspect='auto')
# tickRange = range(0,tmpClust.shape[0])
# xtcks = [x for x in tmpClust.index]
# plt.xticks(tickRange, xtcks,rotation=90)
# plt.yticks(np.arange(len(xtcks)),xtcks)
plt.colorbar()
outF = os.path.join(wkdir, 'dos_clique_heatmap.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

############################################
### Examine cliq connection with members ###
############################################

# positive control (selection of compounds wich are cliq memebers)
n_members = 300
cpSer = pd.Series([item for sublist in cliqFrm['sig'] for item in sublist])
cpSer = cpSer[cpSer.isin(cliqFull.index)]
# iRand = np.random.choice(len(cpLst),n_members,replace=False)
cpRand = np.random.choice(cpSer,n_members,replace=False)
rCliq = cliqFull.reindex(index=cpRand)

# negative control (random selection of compounds)
n_members = 300
summlyCps = anntFrm[anntFrm.pert_type == 'trt_cp'].index.values
cpRand = np.random.choice(summlyCps,n_members,replace=False)
rCliq = cliqFull.reindex(index=cpRand)

# plot distribution
cUnstack = rCliq.unstack()
plt.hist(cUnstack,30)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('median ' + matrixType,fontweight='bold')
plt.title('distribution of DOS-clique connections')
# outF = os.path.join(wkdir, 'positive_control_clique_connection_distribution.png')
outF = os.path.join(wkdir, 'negative_control_clique_connection_distribution.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

### run against dmso and random





### try permuting clique labels --> make different permutations of clique groups

### which DOS compounds have a medain connection greater than 80 lass
maxDos = cliqDos.max(axis=1)
minDos = cliqDos.min(axis=1)
connection_thresh = 85
# positive or negative connection 
# passThresh = (maxDos >= connection_thresh) + (minDos <= -connection_thresh) # above threshold or bellow -threshold
# positive connection only
passThresh = (maxDos >= connection_thresh)
connectedDos = cliqDos.ix[passThresh,:]

### which clique member compounds have a medain connection greater than 85 lass - to a clique
isMember = cliqFull.index.isin(cpSer)
maxCp = cliqFull.max(axis=1)
minCp = cliqFull.min(axis=1)
connection_thresh = 85
# positive or negative connection 
# passThresh = (maxDos >= connection_thresh) + (minDos <= -connection_thresh) # above threshold or bellow -threshold
# positive connection only
passThresh = (maxCp >= connection_thresh)
passThresh = passThresh[isMember]
isConnected = passThresh[passThresh]
connMembers = cliqFull.ix[isConnected.index,:]
# cpSer = cpSer[cpSer.isin(cliqFull.index)]


#########################################
### cluster sig_cliquescore_tool results ###
#########################################

#rows
# mtrx = cliqDos
# mtrx = connectedDos
mtrx = connMembers
Y = scipy.cluster.hierarchy.linkage(mtrx, method='centroid')
Z = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
cOrder = Z['leaves']
iPCL = mtrx.index[cOrder]
clustered = mtrx.reindex(index=iPCL)
#columns 
Y_ = scipy.cluster.hierarchy.linkage(mtrx.T, method='centroid')
Z_ = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
cOrder_ = Z['leaves']
iPCL_columns = mtrx.columns[cOrder]
clustered = mtrx.reindex(index=iPCL,columns=iPCL_columns)

# make heatmap
plt.close()
ccol.set_color_map()
fig = plt.figure(1, figsize=(100, 25))
plt.imshow(clustered.T.values,
    interpolation='nearest',
    aspect='auto')
xtickRange = range(0,clustered.shape[0])
xtcks = [x for x in clustered.index]
ytickRange = range(0,clustered.shape[1])
ytcks = [x for x in clustered.columns]
# cFrm2 = cliqFrm.reindex(clustered.columns)
# ytcks = [x[1]['desc'] + ' - ' + str(x[1]['group_size']) for x in cFrm2.iterrows()]
plt.xticks(xtickRange, xtcks,rotation=90)
plt.yticks(ytickRange, ytcks)
plt.xlabel('compounds')
# plt.yticks(np.arange(len(xtcks)),xtcks)
# plt.title('median connection of DOS to compound class')
plt.title('median connection of clique members to cliques')
plt.colorbar()
# outF = os.path.join(wkdir, 'dos_clique_heatmap_clust.png')
outF = os.path.join(wkdir, 'clique_member_heatmap_thresholded.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

### get pert_inames for cpSer
mc = mu.MongoContainer()
pertInfo = mc.pert_info.find({'pert_id':{'$in':list(cpSer)}},
            {'pert_id':True,'pert_iname':True},toDataFrame=True)



### group x group cluster
groupCorr = np.corrcoef(clustered,rowvar=0)
# make heatmap
plt.close()
fig = plt.figure(1, figsize=(30, 25))
plt.imshow(groupCorr,
    interpolation='nearest',
    aspect='auto',
    cmap=cm.RdBu_r)
ytickRange = range(0,clustered.shape[1])
ytcks = [x for x in clustered.columns]
plt.xticks(ytickRange, ytcks,rotation=90)
plt.yticks(ytickRange, ytcks)
plt.xlabel('compound groups')
plt.title('group-group connection correlation')
plt.colorbar()
outF = os.path.join(wkdir, 'group_by_group_heatmap.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()
