#! /usr/bin/env python

'''
Load in kinase clustering across cell lines

LH 02/2014
'''
import os
import numpy as np, pandas as pd, scipy as sp
from matplotlib import pyplot as plt
from cmap.analytics.statsig import ConnectivitySignificance
from cmap.io import gct
from cmap.util.mongo_utils import CredentialsObject
from pymongo import MongoClient

wkdir = '/xchip/cogs/projects/pharm_class/lhwork/kinase_clustering'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

sFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/kinase_inhibitors/summly/self_rankpt_n302x302.gctx'
sGCT = gct.GCT()
sGCT.read(sFile)
summFrm = sGCT.frame
pInames = sGCT.get_column_meta('pert_iname')
pIDs = sGCT.get_column_meta('pert_id')

sTree = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/kinase_inhibitors/summly/kinase_inhibitor_summly.tree'
# Munge to load in cluster asignments
treeFrm = pd.read_csv(sTree,sep='(')
treeSer = treeFrm.ix[:,0]
treeSplit = treeSer.str.split(':')
treeSplit = treeSplit[~treeSplit.isnull()]
for x in treeSplit:
    print x[0]
listSpace = [x[0] for x in treeSplit]
gSer= pd.Series(listSpace)
inamesCluster = gSer[gSer.isin(pInames)]

#reindex summly results by iname
summIn = summFrm.copy()
summIn.index = pInames
summIn.columns = pInames
summClust = summIn.reindex(index=inamesCluster,columns=inamesCluster)


# load in expected molecular targets of each drug
#import 
aFile = '/xchip/cogs/projects/pharm_class/lhwork/kinase_clustering/drug_annotations.txt'
annFrm = pd.read_csv(aFile,sep='\t')
annFrm = annFrm.reindex(columns=['sum_id','pert_iname','targets'])
annFrm.index = annFrm['pert_iname']
targetFrm = annFrm[annFrm.pert_iname.isin(pInames) & annFrm.sum_id.isin(pIDs)]
targetFrm = targetFrm.reindex(summClust.index)
out = wkdir + '/kinase_targets_clustered.txt'
targetFrm.to_csv(out,sep='\t',headers=True)
plt.close()

# LY-294002
#make heatmap of clustered data
tmpClust = summClust.ix[:30,:30]
fig = plt.figure(1, figsize=(10, 10))
plt.imshow(tmpClust.values,
    interpolation='nearest',
    aspect='auto')
    # cmap=cm.gray_r)
    # vmin=0, 
    # vmax=1,
tickRange = range(0,tmpClust.shape[0])
xtcks = [x for x in tmpClust.index]
plt.xticks(tickRange, xtcks,rotation=90)
plt.yticks(np.arange(len(xtcks)),xtcks)
plt.colorbar()
# plt.xlabel(matrixType + ' threshold')
# plt.ylabel('unique perturbations')
# plt.title('summly false positive rate - based on DMSO')
out = wkdir + '/kinase_summly_clustering.png'
plt.savefig(out, bbox_inches='tight')
plt.close()

# order targets by freq
# loop through each target
# make graph of clustering for each target















