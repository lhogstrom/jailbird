#! /usr/bin/env python

'''
top-level script to compute significance using the statsig module from
cmap.analytics
'''

import numpy as np, pandas as pd, scipy as sp
from matplotlib import pyplot as plt
from cmap.analytics.statsig import ConnectivitySignificance
from cmap.io import gct

# run the DMSO version
score_path = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_n7147x7147.gctx'
# score_path = '/xchip/cogs/projects/connectivity/summly/matrices/asym/matched_lass_n7147x7147.gctx'

null_path = '/xchip/cogs/projects/connectivity/null/dmso/lass_n1000x7147.gctx'
# null_path = '/xchip/cogs/projects/connectivity/null/random/lass_n1000x7147.gctx'

pval_method = 'byrow'
# qval_method = 'byrow'
qval_method = 'byecdf'

# out = '/xchip/cogs/projects/connectivity/null/results_dmso'
# out = '/xchip/cogs/projects/connectivity/null/results_random'
# out = '/xchip/cogs/projects/connectivity/null/results_dmso_assym'
out = '/xchip/cogs/projects/connectivity/null/lh_lass_dmso'

self = ConnectivitySignificance(score_path, null_path, out, 
                                pval_method, qval_method, verbose = True)

self.compute_pvalues()
self.compute_qvalues()
self.generate_outputs()

### q-value by ecdf
reload(ConnectivitySignificance)
self.score_type = 'rnkpt_matched_lass'
self.compute_qvalues_by_ecdf()
# self.compute_qvalues()



### examine Dave's pval calculation:

idx = 'BRD-K70792160'
target = self.data['score'].frame.loc[idx]
null = self.data['null'].frame.loc[idx]
ecdf = ECDF(null)
arg1 = ecdf(target)
arg2 = 1 - ecdf(target)
pvals = 2 * np.minimum(arg1, arg2)
if isinstance(target, pd.Series):
    pvals = pd.Series(pvals, index = target.index, name = 'pvalues')

df = pd.DataFrame(target)
df['ecdf_target'] = arg1
df['inv_target'] = arg2
# df['ecdf_null'] = arg1
df = df.sort(idx)

wkdir = out
fig = plt.figure(1, figsize=(10, 10))
plt.subplot(2,1,1)
a1 = plt.plot(df[idx],df['ecdf_target'],color='b',label='null n=' + str(len(arg1)))
a1 = plt.plot(df[idx],df['inv_target'],color='g',label='null inverse n=' + str(len(arg1)))
# a3 = plt.plot(vals,dEval,color='r',label='DMSO n=' + str(len(dmsoVec))) #
plt.legend(loc=2)
plt.ylabel('F(x)',fontweight='bold')
# plt.xlabel(matrixType,fontweight='bold')
plt.title('ecdf for summly row - ' + idx)
plt.subplot(2,1,2)
h1 = plt.hist(target,30,color='b',range=[-100,100],label=['observed'],alpha=.4,normed=True)
h2 = plt.hist(null,30,color='r',range=[-100,100],label='DMSO',alpha=.3,normed=True)
# plt.legend()
plt.ylabel('freq',fontweight='bold')
plt.xlabel(matrixType,fontweight='bold')
outF = os.path.join(wkdir, idx + '_ecdf.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

