#! /usr/bin/env python

'''
access the signficance of the observed clique results - compared to DMSO
'''

import os
import numpy as np, pandas as pd, scipy as sp
from matplotlib import pyplot as plt
import copy
from matplotlib import cm
from statsmodels.distributions import ECDF

from cmap.analytics.statsig import ConnectivitySignificance
from cmap.io import gct
from cmap.io import gmt
import cmap.util.progress as update

# wkdir = '/xchip/cogs/projects/connectivity/null/clique_analysis/shared_targets_vs_dmso_null'
wkdir = '/xchip/cogs/projects/connectivity/null/clique_analysis/meta_compound_classes_vs_dmso_null'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

# Load Steven's cliques
grpMin = 2
# rFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/cpd_targets_n268/summly/self_rankpt_n342x342.gctx'
# cFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/cpd_targets_n268/summly/signature_info.txt'
rFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/cpd_groups_n147/summly/self_rankpt_n1096x1096.gctx'
cFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/cpd_groups_n147/summly/signature_info.txt'

classFrm = pd.read_csv(cFile,sep='\t')
classGrp = classFrm.groupby('group_id')
grpDict = {}
for grp in classGrp.groups:
    igrp = classGrp.groups[grp]
    grpFrm = classFrm.reindex(igrp)
    pIds = list(grpFrm['pert_id'])
    if len(pIds) < grpMin:
        continue
    grpDict[grp] = pIds
grpSer = pd.Series(grpDict)
grpSer.name = 'sig'
cliqueLabels = pd.DataFrame(grpSer)
cliqueLabels['id'] = cliqueLabels.index

cList = [item for sublist in cliqueLabels['sig'] for item in sublist]
cSet = set(cList)

# load observed score data
# thresholded
# rFile = '/xchip/cogs/projects/connectivity/null/clique_analysis/dmso_q_thresholded_asym_lass_matrix/jan28/my_analysis.sig_cliqueselect_tool.2014012814320559/summly/self_rankpt_n379x379.gctx'
# non-thresholded asym
# rFile = '/xchip/cogs/projects/connectivity/null/clique_analysis/baseline_lass_asym_matrix/jan28/my_analysis.sig_cliqueselect_tool.2014012814364180/summly/self_rankpt_n379x379.gctx'
gt1 = gct.GCT()
gt1.read(rFile)
sFrm = gt1.frame
sFrm.columns = gt1.get_column_meta('pert_id')
#check that all clique members are in the observed matrix
if not (sFrm.index.isin(cSet)).all():
    print "not all clique data loaded"

# load null 
dFile = '/xchip/cogs/projects/connectivity/null/dmso/lass_n1000x7147.gctx'
gt = gct.GCT(dFile)
gt.read()
dmsoFrm = gt.frame
dmsoFrm.columns = gt.get_column_meta('id')
dmsoCM = dmsoFrm[dmsoFrm.index.isin(cSet)]
rowMedian = dmsoCM.median(axis=1)

def no_diagonal_unstack(frm):
    'return an unstacked matrix without the diagonal'
    np.fill_diagonal(frm.values, np.nan)
    overlapSer = frm.unstack()
    overlapSer = overlapSer[~overlapSer.isnull()] #remove nulls 
    return overlapSer

### compare observed to null
#construct 
graph=True
pvalDict = {}
progress_bar = update.DeterminateProgressBar('group p-val computation')
for iicliq,icliq in enumerate(cliqueLabels.index):
    progress_bar.update('count', iicliq, len(cliqueLabels.index))
    cName = cliqueLabels.ix[icliq,'id']
    pIds = cliqueLabels.ix[icliq,'sig']
    smFrm = sFrm.reindex(index=pIds,columns=pIds)
    uFrm = no_diagonal_unstack(smFrm)
    medObs = uFrm.median()
    rMed = rowMedian[pIds]
    fig = plt.figure(1, figsize=(10, 10))
    # make matrix of equal size using null
    nperm = 10000
    permDict = {}
    for iperm in range(nperm):
        iRand = np.random.choice(range(0,dmsoFrm.shape[1]),size=(len(pIds)))
        iRandCol = dmsoFrm.columns[iRand] #random column names
        smDmso = dmsoFrm.reindex(index=pIds,columns=iRandCol)
        # remove identity cells and unstack
        uDmso = no_diagonal_unstack(smDmso)
        medDmso = uDmso.median()
        permDict[iperm] = medDmso
    nullSer = pd.Series(permDict)
    #two tailed p-value
    ecdf = ECDF(nullSer)
    arg1 = ecdf(medObs)
    arg2 = 1 - ecdf(medObs)
    pval = 2 * np.minimum(arg1, arg2)
    #set p-val min
    if pval == 0:
        pval = 1/float(nperm)
    pvalDict[cName] = pval
    if graph:
        # graph heatmap of each
        plt.imshow(smFrm.values,
            interpolation='nearest',
            aspect='auto',
            vmin=-100, 
            vmax=100,
            cmap=cm.RdBu_r)
        tickRange = range(0,smFrm.shape[1])
        xtcks = [x for x in smFrm.columns]
        plt.xticks(tickRange, xtcks)
        plt.yticks(tickRange,xtcks)
        plt.colorbar()
        plt.title(cName)
        out = wkdir + '/' + cName + '_observed.png'
        plt.savefig(out, bbox_inches='tight')
        plt.close()
        # median null hist
        h2 = plt.hist(nullSer,30,color='b',range=[-100,100],label='DMSO',alpha=.8)
        h1 = plt.hist(np.array([medObs]),30,color='r',range=[-100,100],label=['observed'],alpha=.3,weights=[1000])
        plt.legend()
        plt.xlabel('lass scores',fontweight='bold')
        plt.title('clique median score - ' + cName)
        out = wkdir + '/' + cName + '_null_medians.png'
        plt.savefig(out, bbox_inches='tight')
        plt.close()
pvalSer = pd.Series(pvalDict)
pvalSer.name = 'group_median_pval'
outF = wkdir + '/clique_group_median_pval.txt'
pvalSer.to_csv(outF, sep='\t',header=True)

#plot row median of DMSO
plt.hist(rowMedian,30)
plt.xlabel('lass scores',fontweight='bold')
plt.ylabel('count',fontweight='bold')
plt.title('clique members - dmso row median')
outF = wkdir + '/dmso_row_median.png'
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

