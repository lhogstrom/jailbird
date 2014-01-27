#! /usr/bin/env python

'''
-threshold lass matrix acording to connection q-val
-prep to run Rajiv's sig_cliqueselect_tool in matlab
'''

import numpy as np, pandas as pd, scipy as sp
from matplotlib import pyplot as plt
from cmap.analytics.statsig import ConnectivitySignificance
from cmap.io import gct

# matched lass connection matrix
lFile = '/xchip/cogs/projects/connectivity/null/clique_analysis/baseline_lass_matrix/matched_lass_n7147x7147.gctx'
gt = gct.GCT()
gt.read(lFile)
lassMtrx = gt.frame

# q-value matrix from fdr correction
qFile = '/xchip/cogs/projects/connectivity/null/results_dmso_sym/qvalues_n7147x7147.gctx'
# qFile = '/xchip/cogs/projects/connectivity/null/results_random_sym/qvalues_n7147x7147.gctx'
gt2 = gct.GCT()
gt2.read(qFile)
qMtrx = gt2.frame

if not (lassMtrx.index == qMtrx.index).all():
    print 'matrix indices do not match'
if not (lassMtrx.columns == qMtrx.columns).all():
    print 'matrix indices do not match'

#chech that index and columns of matrices are the same

qMask = qMtrx.copy()
qThresh = .1
islt = qMtrx <= qThresh
qMask[islt] = 1
qMask[~islt] = 0

lassMasked = qMask*lassMtrx
ogt = gct.GCT()
ogt.build_from_DataFrame(lassMasked)
outF = '/xchip/cogs/projects/connectivity/null/clique_analysis/q_thresholded_lass_matrix/matched_lass_n7147x7147.gctx'
ogt.write(outF)

#load in thresholded matrix
gt3 = gct.GCT()
gt3.read(outF)
outFrm = gt3.frame

# run matlab code like this:
# sig_cliqueselect_tool('clique','/xchip/cogs/projects/pharm_class/cp_cliques_current.gmt', 'inpath', '/xchip/cogs/projects/connectivity/summly/matrices/','out','/xchip/cogs/projects/connectivity/summly/matrices/')