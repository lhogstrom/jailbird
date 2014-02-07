#! /usr/bin/env python

'''
Load in kinase clustering across cell lines

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
targetFrm.targets = targetFrm.targets.str.split(', ')
targetFrm = targetFrm.reindex(summClust.index)
out = wkdir + '/kinase_targets_clustered.txt'
targetFrm.to_csv(out,sep='\t',headers=True)

# order targets by freq
# loop through each target
tarLst = [item for sublist in targetFrm.targets for item in sublist]
tarSer = pd.Series(tarLst)
# targetSet = set(tarLst)
tarCounts = tarSer.value_counts()
targetMatch = targetFrm.copy()
nTarget = targetMatch.shape[0]
boolFrm = pd.DataFrame(np.zeros([nTarget,len(tarCounts)]),index=targetFrm.index,columns=tarCounts.index)
for tar in tarCounts.index:
    def check_cell_contents(x):
        if tar in x:
            return 1
        else:
            return 0
    boolMatch = targetFrm.targets.apply(check_cell_contents)
    targetMatch[tar] = boolMatch
    boolFrm.ix[:,tar] = boolMatch
out = wkdir + '/kinase_targets_match.txt'
targetMatch.to_csv(out,sep='\t',headers=True)

#make heatmap of clustered data
tmpClust = summClust.ix[:30,:30]
fig = plt.figure(1, figsize=(10, 10))
plt.subplot(1,2,1)
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
# target status
plt.subplot(1,2,2)
smFrm = pd.DataFrame(boolFrm['MAP2K1'][:30])
plt.imshow(smFrm,
    interpolation='nearest',
    aspect='auto')
#save
out = wkdir + '/kinase_summly_clustering.png'
plt.savefig(out, bbox_inches='tight')
plt.close()

### target status heatmap
#make heatmap of clustered data
nStart = 40
nFinish = 65
tmpClust = summClust.ix[nStart:nFinish,nStart:nFinish]
fig = plt.figure(1, figsize=(10, 10))
gs = gridspec.GridSpec(1, 2, width_ratios=[40, 1]) 
ax0 = plt.subplot(gs[0])
# plt.subplot(1,2,1)
ax0.imshow(tmpClust.values,
    interpolation='nearest',
    aspect='auto',
    vmin=-100, 
    vmax=100)
tickRange = range(0,tmpClust.shape[0])
xtcks = [x for x in tmpClust.index]
plt.xticks(tickRange, xtcks,rotation=90)
plt.yticks(np.arange(len(xtcks)),xtcks)
# plt.colorbar()
# plt.xlabel(matrixType + ' threshold')
# plt.ylabel('unique perturbations')
plt.title('Kinase clustering')
# target status
target = 'AURKB'
smFrm = pd.DataFrame(boolFrm[target][nStart:nFinish])
ax1 = plt.subplot(gs[1])
ax1.imshow(smFrm,
    interpolation='nearest',
    aspect='auto')
ax1.axes.get_yaxis().set_visible(False)
plt.xticks([0], ['labeled target - ' + target],rotation=90)
out = wkdir + '/' + target + '_target_label_clustering.png'
plt.savefig(out, bbox_inches='tight')
plt.close()

### rolling sum window
for tar in tarCounts.index[:70]:
    target = tar
    window=10
    rollSum = pd.stats.moments.rolling_sum(boolFrm[target],window)
    out = wkdir + '/' + target + '_rolling_sum.png'
    plt.plot(rollSum)
    plt.ylim((0,window))
    plt.xlabel('cluster axis')
    plt.ylabel('rolling sum')   
    plt.title(target + ' - target density from clustering - window = ' +str(window)) 
    plt.savefig(out, bbox_inches='tight')
    plt.close()


# heatmap of targets
tmpBool = boolFrm.ix[:,:30]
fig = plt.figure(1, figsize=(10, 50))
plt.imshow(tmpBool.values,
    interpolation='nearest',
    aspect='auto',
    cmap=cm.Greys)
    # vmin=0, 
    # vmax=1,
ytickRange = range(0,tmpBool.shape[0])
xtickRange = range(0,tmpBool.shape[1])
ytcks = [x for x in tmpBool.index]
xtcks = [x for x in tmpBool.columns]
# ax.xaxis.set_ticks_position('top')
plt.xticks(xtickRange, xtcks,rotation=90)
plt.yticks(ytickRange,ytcks)
# plt.colorbar()
# plt.xlabel(matrixType + ' threshold')
# plt.ylabel('unique perturbations')
# plt.title('summly false positive rate - based on DMSO')
out = wkdir + '/kinase_bool_targets.png'
plt.savefig(out, bbox_inches='tight')
plt.close()

# loop through each cell line - apply summly dendrogram to order to all
# make heatmap




