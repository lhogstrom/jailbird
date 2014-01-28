#! /usr/bin/env python

'''
-threshold lass matrix acording to connection q-val
-prep to run Rajiv's sig_cliqueselect_tool in matlab
'''

import os
import numpy as np, pandas as pd, scipy as sp
from matplotlib import pyplot as plt
from cmap.analytics.statsig import ConnectivitySignificance
from cmap.io import gct
import copy

wkdir = '/xchip/cogs/projects/connectivity/null/clique_analysis/random_q_thresholded_asym_lass_matrix'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

# matched lass connection matrix
# lFile = '/xchip/cogs/projects/connectivity/null/clique_analysis/baseline_lass_sym_matrix/matched_lass_n7147x7147.gctx'
lFile = '/xchip/cogs/projects/connectivity/null/clique_analysis/baseline_lass_asym_matrix/matched_lass_n7147x7147.gctx'
gt = gct.GCT()
gt.read(lFile)
lassMtrx = gt.frame

# q-value matrix from fdr correction
# qFile = '/xchip/cogs/projects/connectivity/null/results_dmso_sym/qvalues_n7147x7147.gctx'
# qFile = '/xchip/cogs/projects/connectivity/null/results_random_sym/qvalues_n7147x7147.gctx'
# qFile = '/xchip/cogs/projects/connectivity/null/results_dmso_assym/qvalues_n7147x7147.gctx'
qFile = '/xchip/cogs/projects/connectivity/null/results_random_assym/qvalues_n7147x7147.gctx'
gt2 = gct.GCT()
gt2.read(qFile)
qMtrx = gt2.frame

if not (lassMtrx.index == qMtrx.index).all():
    print 'matrix indices do not match'
if not (lassMtrx.columns == qMtrx.columns).all():
    print 'matrix indices do not match'

#chech that index and columns of matrices are the same
qThresh = .2
isgt = qMtrx > qThresh
lassMasked = lassMtrx.copy()
lassMasked[isgt] = 0

ogt = gct.GCT()
gt.mk_cdesc()
gt.mk_rdesc()
ogt.build_from_DataFrame(lassMasked,rdesc=gt.rdesc,cdesc=gt.cdesc)
# ogt = copy.copy(gt)
# ogt.matrix = lassMasked.values
#replace matrix
outF = wkdir + '/matched_lass_n7147x7147.gctx'
ogt.write(outF)

#load in thresholded matrix
gt3 = gct.GCT()
gt3.read(outF)
outFrm = gt3.frame

#plot distribution of significant lass scores:
g = lassMasked[~(lassMasked == 0)]
scoreList = g.unstack()
scoreList = scoreList[~np.isnan(scoreList)]

plt.hist(scoreList,30)
plt.xlabel('lass scores',fontweight='bold')
plt.ylabel('count',fontweight='bold')
plt.title('lass random significance - q<.2')
outF = wkdir + '/random_significant_lass_scores_hist.png'
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

# what does baseline look like?
lassList = lassMtrx.unstack()
lassList = lassList[~np.isnan(lassList)]
plt.hist(lassList,30)
plt.xlabel('lass scores',fontweight='bold')
plt.ylabel('count',fontweight='bold')
plt.title('all lass scores')
outF = wkdir + '/all_lass_scores_hist.png'
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()


# what do the random matrices look like?

# run matlab code like this:
# sig_cliqueselect_tool('clique','/xchip/cogs/projects/pharm_class/cp_cliques_current.gmt', 'inpath', '/xchip/cogs/projects/connectivity/summly/matrices/','out','/xchip/cogs/projects/connectivity/summly/matrices/')

# clique by null --> for what is the median score for a column based on a clique
