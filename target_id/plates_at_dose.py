 #! /usr/bin/env python
'''
examine/search for plates with dose data
'''
import glob, HTML
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection as fdr
import cmap.util.mongo_utils as mu
from cmap.tools import sig_slice_tool
from cmap.io import gct,plategrp,rnk
import cmap.util.progress as progress
import subprocess
import datetime
import cmap.util.tool_ops as to
import cmap.analytics.dgo as dgo

#search through plates with dose data
CM = mu.CMapMongo()
hogList = CM.find({'sig_id':{'$regex':'HOG'}},{'pert_id':True})
hogSet = set(hogList)

CM = mu.CMapMongo()
pclbList = CM.find({'sig_id':{'$regex':'PCLB'}},{'pert_id':True})
pclbSet = set(pclbList)

CM = mu.CMapMongo()
brafList = CM.find({'sig_id':{'$regex':'BRAF'},'pert_type':'trt_cp'},{'pert_id':True})
brafSet = set(brafList)

# lSet = [hogSet,pclbSet,brafSet]
lSet = [pclbSet,brafSet]
combinedSet =frozenset().union(*lSet)

#check how many doses were recorded for each cp
pDescDict = {}
doseDict = {}
for brd in pclbSet:
    iQ = CM.find({'sig_id':{'$regex':'PCLB'},'pert_id':brd},{'pert_dose':True,'pert_iname':True})
    # or expression - differnt plates
    # iQ = CM.find({'pert_id':brd, '$or':[{'sig_id':{'$regex':'PCLB'}},{'sig_id':{'$regex':'BRAF'}}]},{'pert_dose':True,'pert_iname':True})
    idoses = [x['pert_dose'] for x in iQ]
    doseDict[brd] = set(idoses)
    inames = [x['pert_iname'] for x in iQ]
    pDescDict[brd] = set(inames)
hogSet.remove('CMAP-HSF-HOGA1')
hogSet.remove('DMSO')

CM.find({'$or':[{'pert_id':p1},{'pert_id':p2}]},{'pert_dose':True,'pert_iname':True})



