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

wkdir = '/xchip/cogs/projects/pharm_class/lhwork/summly_clustering_16Feb2014'
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

sTree = '/xchip/cogs/projects/connectivity/clustering/matched_lass_dendro.tre'
# Munge to load in cluster asignments
treeFrm = pd.read_csv(sTree,sep='(')
treeSer = treeFrm.ix[:,0]
treeSplit = treeSer.str.split(':')
treeSplit = treeSplit[~treeSplit.isnull()]
listSpace = [x[0] for x in treeSplit]
gSer= pd.Series(listSpace)
inamesCluster = gSer[gSer.isin(pInameType)]
clustFrm = anntFrm.reindex(inamesCluster.values)
clustFrm['order'] = range(clustFrm.shape[0])

###############################
### load groupings         #### 
###############################

cFile = '/xchip/cogs/projects/pharm_class/rnwork/cliques/cpd_groups_n147.gmt'
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
window=10
prog = progress.DeterminateProgressBar('cliq graph')
for icliq,cliq in enumerate(cliqFrm.desc):
    prog.update(cliq,icliq,len(cliqFrm.desc))
    cliqMod = cliqFrm.ix[icliq,'desc_mod']
    brds = cliqFrm.ix[icliq,'sig']
    boolSer = clustFrm.pert_id.isin(brds)
    rollSum = pd.stats.moments.rolling_sum(boolSer,window)
    out = wkdir + '/' + cliqMod + '_rolling_sum.png'
    plt.plot(rollSum)
    plt.ylim((0,window))
    plt.xlabel('cluster axis')
    plt.ylabel('rolling sum')   
    plt.title(cliq + ' - target density from clustering - window = ' +str(window)) 
    plt.savefig(out, bbox_inches='tight')
    plt.close()
