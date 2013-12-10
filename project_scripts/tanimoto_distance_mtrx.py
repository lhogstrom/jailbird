'''
-use tanimoto matrix
-see if there is a relationship between structure and expression
-1) across the whole matrix
-2) in pharmalogical classes
'''
import cmap.analytics.pcla as pcla
import numpy as np
import os
import cmap
import pandas as pd
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu

# load in summly matrix
rnkptGCT = gct.GCT()
rnkptGCT.read('/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_sym_n7322x7322.gctx')
rnkpt = rnkptGCT.frame
rnkpt.columns = rnkpt.index
summBrds = rnkpt.index.values

# load in tanimoto matrix
tFile = '/xchip/cogs/jgould/tanimoto.gctx'
gt = gct.GCT()
gt.read(tFile)
# gt.read(src=tFile,cid=list(summBrds),rid=list(summBrds))
tanFrm = gt.frame

tanMtch = tanFrm.reindex(index=rnkpt.index,columns=rnkpt.index)

# make array of tanimoto distances and summly scores
nm = rnkpt.shape[0]
iUp = np.tril_indices(nm)

iUp = np.tril_indices(100)
rnkptCopy = rnkpt.copy().ix[:100,:100]
rnkptCopy.values[iUp] = np.nan
rnkptArray = rnkptCopy.unstack().values
rnkptArray = rnkptArray[~np.isnan(rnkptArray)]

tanCopy = tanMtch.copy().ix[:100,:100]
tanCopy.values[iUp] = np.nan
tanArray = tanCopy.unstack().values
tanArray = tanArray[~np.isnan(tanArray)]


#how many nans are in tanimoto matrix?
tanCopy = tanMtch.copy()
tanArray = tanCopy.unstack().values
sumNan = sum(np.isnan(tanArray))


