'''
-examine tanimoto matrix
-compare to summly matrix

Larson Hogstrom, 4/2014
'''
import cmap.io.gct as gct
import pandas as pd
import matplotlib.pyplot as plt
import os
import cmap.io.gmt as gmt
import cmap.util.progress as update
import cmap.util.mongo_utils as mu

Ftanimoto = '/xchip/cogs/jgould/tanimoto.gctx'
gt = gct.GCT()
gt.read(Ftanimoto)
tanimoto = gt.frame

g = tanimoto.ix[:,'BRD-K81418486']
g = g[~g.isnull()]
plt.hist(g,30)