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

wkdir = '/xchip/cogs/projects/connectivity/null/clique_analysis/clique_vs_dmso_null'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

# load cliques
classGMT = '/xchip/cogs/projects/pharm_class/cp_cliques_current.gmt'
gmtDict = gmt.read(classGMT)
cliqueLabels = pd.DataFrame(gmtDict)
# create set of all clique members
cList = [item for sublist in cliqueLabels['sig'] for item in sublist]
cSet = set(cList)

# load observed score data
# thresholded
# rFile = '/xchip/cogs/projects/connectivity/null/clique_analysis/dmso_q_thresholded_asym_lass_matrix/jan28/my_analysis.sig_cliqueselect_tool.2014012814320559/summly/self_rankpt_n379x379.gctx'
# non-thresholded asym
rFile = '/xchip/cogs/projects/connectivity/null/clique_analysis/baseline_lass_asym_matrix/jan28/my_analysis.sig_cliqueselect_tool.2014012814364180/summly/self_rankpt_n379x379.gctx'
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
graph=False
pvalDict = {}
progress_bar = update.DeterminateProgressBar('group p-val computation')
for icliq in cliqueLabels.index:
    progress_bar.update('count', icliq, len(cliqueLabels.index))
    cName = cliqueLabels.ix[icliq,'id']
    pIds = cliqueLabels.ix[icliq,'sig']
    smFrm = sFrm.reindex(index=pIds,columns=pIds)
    uFrm = no_diagonal_unstack(smFrm)
    medObs = uFrm.median()
    rMed = rowMedian[pIds]
    fig = plt.figure(1, figsize=(10, 10))
    # make matrix of equal size using null
    nperm = 1000
    permDict = {}
    for iperm in range(nperm):
        iRand = np.random.choice(range(0,dmsoFrm.shape[1]),size=(len(pIds)))
        iRandCol = dmsoFrm.columns[iRand] #random column names
        smDmso = dmsoFrm.reindex(index=pIds,columns=iRandCol)
        # remove identity cells and unstack
        uDmso = no_diagonal_unstack(smDmso)
        medDmso = uDmso.median()
        permDict[iperm] = medDmso
    medianSer = pd.Series(permDict)
    #two tailed p-value
    ecdf = ECDF(medianSer)
    arg1 = ecdf(medObs)
    arg2 = 1 - ecdf(medObs)
    pval = 2 * np.minimum(arg1, arg2)
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
        # plt.xlabel(matrixType + ' threshold')
        # plt.ylabel('unique perturbations')
        plt.title(cName)
        out = wkdir + '/' + cName + '_observed.png'
        plt.savefig(out, bbox_inches='tight')
        plt.close()
        # median null hist
        plt.hist(medianSer,30)
        plt.hist(medObs)
        plt.title(cName)
        out = wkdir + '/' + cName + '_null_medians.png'
        plt.savefig(out, bbox_inches='tight')
pvalSer = pd.Series(pvalDict)





    





#plot row median of DMSO
plt.hist(rowMedian,30)
plt.xlabel('lass scores',fontweight='bold')
plt.ylabel('count',fontweight='bold')
plt.title('clique members - dmso row median')
outF = wkdir + '/dmso_row_median.png'
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

