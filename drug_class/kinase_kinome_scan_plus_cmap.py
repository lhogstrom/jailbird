#! /usr/bin/env python

'''
load in kinome scan data - compare to CMAP

LH 02/2014
'''

import os
import numpy as np, pandas as pd, scipy as sp
import glob
from matplotlib import pyplot as plt
from matplotlib import cm
from statsmodels.distributions import ECDF

import cmap.analytics.statsig as statsig
from cmap.io import gct
from cmap.io import gmt

wkdir = '/xchip/cogs/projects/pharm_class/lhwork/kinase_clustering'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

## in normalized matrix - negative values mean sensitivity to a given kinase
kFile = '/xchip/cogs/stacks/STK034_HMSKSCAN/STK034_zs_uniq_nonan_n106x355.gctx'
## two types Kd and percent control
# nonNormalized = '/xchip/cogs/stacks/STK034_HMSKSCAN/STK034_n165x515.gctx'
kGCT = gct.GCT()
kGCT.read(kFile)
kFrm = kGCT.frame
pInames = kGCT.get_column_meta('pert_iname')
pIDs = kGCT.get_column_meta('pert_id')
kFrm.columns = pInames

# load clustering tree
ktFile = '/xchip/cogs/stacks/STK034_HMSKSCAN/figures/kscan_n84.tree'
# Munge to load in cluster asignments
treeFrm = pd.read_csv(ktFile,sep='(')
treeSer = treeFrm.ix[:,0]
treeSplit = treeSer.str.split(':')
treeSplit = treeSplit[~treeSplit.isnull()]
listSpace = [x[0] for x in treeSplit]
gSer= pd.Series(listSpace)
inamesCluster = gSer[gSer.isin(pInames)]

# load in expected molecular targets of each drug
aFile = '/xchip/cogs/projects/pharm_class/lhwork/kinase_clustering/drug_annotations.txt'
annFrm = pd.read_csv(aFile,sep='\t')
annFrm = annFrm.reindex(columns=['sum_id','pert_iname','targets'])
annFrm.index = annFrm['pert_iname']
targetFrm = annFrm[annFrm.pert_iname.isin(inamesCluster) & annFrm.sum_id.isin(pIDs)]
targetFrm.targets = targetFrm.targets.str.split(', ')
# targetFrm = targetFrm.reindex(inamesCluster.index)

# re-index
clustFrm = kFrm.reindex(columns=inamesCluster.values)
# testDrug = 'erlotinib'
# sortFrm = clustFrm.sort(testDrug)

#check target annotations against kinome data
rnkDict = {}
rnkList = []
for drug in clustFrm.columns:
    drugSer = kFrm.ix[:,drug]
    drugSer = drugSer.order()
    drugRank = drugSer.rank(method='first')
    targets = targetFrm.ix[drug,'targets']
    targetRank = drugRank.reindex(targets)
    # targetRank = targetRank[~targetRank.isnull()]
    tarRnkPair = zip(targetRank.index, targetRank.values)
    rnkDict[drug] = tarRnkPair
    rnkValues = targetRank[~targetRank.isnull()].values
    rnkList.extend(rnkValues)
rnkSer = pd.Series(rnkDict)
rnkSer.index.name = 'drug'
rnkSer.name = 'target_ranks'
out = wkdir + '/kinome_scan_cmap_target_ranks.txt'
rnkSer.to_csv(out,sep='\t',header=True)

# plot target rankings
plt.hist(rnkList,30)
out = wkdir + '/kinome_scan_cmap_target_ranks.png'
plt.xlabel('kscan rank')
plt.ylabel('freq')   
plt.title('Kinome scan results for expected drug targets, n=84 drugs') 
plt.savefig(out, bbox_inches='tight')
plt.close()














