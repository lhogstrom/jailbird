
#! /usr/bin/env python
'''
analyze the DOSBIO plates - combine with data in cmap database to perform query

use DOS signatures generate queries of the CGS data (cell line specific results)
'''

import os
import cmap.io.gct as gct
import cmap.analytics.sc as sc
import glob as glob
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import cmap.analytics.dgo as dgo

### get DOS BIO cps from mongo
CM = mutil.CMapMongo()
# pert_List = CM.find({'sig_id':{'$regex':'DOSBIO'},'pert_iname':{'$regex':'BRD'}},{'sig_id':True,'pert_id':True,'pert_iname':True})
pert_List = CM.find({'sig_id':{'$regex':'DOSBIO'},'pert_iname':{'$regex':'BRD'}},{'pert_id':True})
dosbioSet = set(pert_List)
# check to make sure the brds are DOS compounds and don't represent known compounds
for brd in dosbioSet:
	





