#! /usr/bin/env python

'''
access the signficance of the observed clique results - compared to DMSO
'''

import os
import numpy as np, pandas as pd, scipy as sp
from matplotlib import pyplot as plt
from cmap.analytics.statsig import ConnectivitySignificance
from cmap.io import gct
from cmap.io import gmt
import copy



# load cliques
classGMT = '/xchip/cogs/projects/pharm_class/cp_cliques_current.gmt'
gmtDict = gmt.read(classGMT)
cliqueLabels = pd.DataFrame(gmtDict)
# create set of all clique members
cList = [item for sublist in cliqueLabels['sig'] for item in sublist]
cSet = set(cList)

# load observed score data
rFile = '/xchip/cogs/projects/connectivity/null/clique_analysis/dmso_q_thresholded_asym_lass_matrix/jan28/my_analysis.sig_cliqueselect_tool.2014012814320559/summly/self_rankpt_n379x379.gctx'
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



