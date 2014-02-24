#! /usr/bin/env python

'''
-Load summly clustering
-Evaluate 

LH 02/2014
'''
import os
import numpy as np, pandas as pd, scipy as sp
from matplotlib import pyplot as plt
from pymongo import MongoClient
from cmap.analytics.statsig import ConnectivitySignificance
from cmap.io import gct
from cmap.util.mongo_utils import CredentialsObject
from matplotlib import cm
import matplotlib.gridspec as gridspec
import cmap.io.gmt as gmt
import cmap.util.progress as progress
import cmap.util.mongo_utils as mu

wkdir = '/xchip/cogs/projects/pharm_class/lhwork/summly_clustering_24Feb2014'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

###########################
### load summly matrix ####
###########################

### use lass matched matrix
sFile = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_n7147x7147.gctx'
# sFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/kinase_inhibitors/summly/self_rankpt_n302x302.gctx'
sGCT = gct.GCT()
sGCT.read(sFile)
summFrm = sGCT.frame
pInames = sGCT.get_column_meta('pert_iname')
pIDs = sGCT.get_column_meta('pert_id')
pType = sGCT.get_column_meta('pert_type')
pInameType = []
for i,x in enumerate(pInames):
    pInameType.append(pInames[i]+ '.' +pType[i])
anntFrm = pd.DataFrame({'pert_id':pIDs,'pert_type':pType,'pert_iname':pInames},index=pInameType)

###############################
### load hierarchical tree #### 
###############################

### load XML file 
sTreeXML = '/xchip/cogs/projects/connectivity/clustering/matched_lass_dendro.xml'
tstPanLst = open(sTreeXML).read().splitlines()
labelList = []
for str1 in tstPanLst:
    if 'label="' in str1:
        s1 = str1.split('label="')[1]
        s2 = s1.split('">')[0]
        labelList.append(s2)
labelList.pop(0)
clustFrm = anntFrm.reindex(labelList)
clustFrm['order'] = range(clustFrm.shape[0])
outF = wkdir + '/summly_dendrogram_table.txt'
clustFrm.to_csv(outF,sep='\t')

##########################################
### load DOS compounds in summly space ###
##########################################
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

dosBrds = get_dos_BRDs()
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_n7147x7147.gctx'
matrixType = 'rnkpt_matched_lass'
dosSer = get_summly_dos_indeces(dosBrds,mtrxSummly)

###############################
### load groupings         #### 
###############################

# cFile = '/xchip/cogs/projects/pharm_class/rnwork/cliques/cpd_groups_n147.gmt'
cFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140213/cliques.gmt'
cliqueGMT = gmt.read(cFile)
cliqFrm = pd.DataFrame(cliqueGMT)
cliqFrm['desc_mod'] = cliqFrm['desc']
cliqFrm['desc_mod'] = cliqFrm['desc_mod'].str.replace("/","-")
cliqFrm['desc_mod'] = cliqFrm['desc_mod'].str.replace(" ","_")
cliqFrm['desc_mod'] = cliqFrm['desc_mod'].str.replace("&","_")
cliqFrm['desc_mod'] = cliqFrm['desc_mod'].str.replace("?","_")
cliqFrm['desc_mod'] = cliqFrm['desc_mod'].str.replace("(","_")
cliqFrm['desc_mod'] = cliqFrm['desc_mod'].str.replace(")","_")
cliqFrm['desc_mod'] = cliqFrm['desc_mod'].str.replace("'","_")
cliqFrm['desc_mod'] = cliqFrm['desc_mod'].str.replace('\xce','_')
cliqFrm['desc_mod'] = cliqFrm['desc_mod'].str.replace('\xba','_')
cliqFrm['desc'] = cliqFrm['desc'].str.replace('\xce','_')
cliqFrm['desc'] = cliqFrm['desc'].str.replace('\xba','_')
### rolling sum window
graphDir = wkdir + '/window_min'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
window=10
group_min = 4 # minimum clique members within a window
clustered_groups = {}
dosDict = {}
prog = progress.DeterminateProgressBar('cliq graph')
for icliq,cliq in enumerate(cliqFrm.desc):
    prog.update(cliq,icliq,len(cliqFrm.desc))
    cliqMod = cliqFrm.ix[icliq,'desc_mod']
    brds = cliqFrm.ix[icliq,'sig']
    boolSer = clustFrm.pert_id.isin(brds)
    isDos = clustFrm.pert_id.isin(dosSer.index)
    rollSum = pd.stats.moments.rolling_sum(boolSer,window)
    rollSum.name = 'rolling_group_count'
    nrollSum = rollSum[~rollSum.isnull()]
    rollFrm = pd.DataFrame(rollSum)
    rollFrm['location'] = np.arange(len(rollSum))
    rollFrm[cliqMod] = boolSer
    rollFrm['is_dos'] = isDos
    rollFrm['pert_id'] = clustFrm.pert_id
    if max(nrollSum) > group_min:
        peakSer = nrollSum[nrollSum > group_min]
        clustered_groups[cliq] = peakSer
        rollFrm.ix[nrollSum[nrollSum > group_min].index,:]
        rMin = rollFrm.index.get_loc(peakSer.index[0])
        rMax = rollFrm.index.get_loc(peakSer.index[-1])
        rollRange = np.arange(rMin,rMax)
        if len(rollRange) <= 100: # skip if cluster is too long
            localFrm = rollFrm.ix[(rMin-10):(rMax+10),:]
            out = graphDir + '/' + cliqMod + '_local_dendrogram_table.txt'
            localFrm.to_csv(out,sep='\t')
            # any local dos compounds
            if localFrm.is_dos.any():
                dosDict[cliq] = list(localFrm[localFrm.is_dos].pert_id)
        # create graph
        out = graphDir + '/' + cliqMod + '_rolling_sum.png'
        plt.plot(rollSum)
        plt.ylim((0,window))
        plt.xlabel('cluster axis')
        plt.ylabel('rolling sum')   
        plt.title(cliq + ' - target density from clustering - window = ' +str(window)) 
        plt.savefig(out, bbox_inches='tight')
        plt.close()

# list of DOS compounds that are close by clique region of the dendrogram
dosDendro = pd.Series(dosDict)
outF = wkdir + '/DOS_cps_proximal_to_local_cliques.txt'
dosDendro.to_csv(outF,sep='\t')


# dosConnected = [item for sublist in dosDict.values() for item in sublist]



