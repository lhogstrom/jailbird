'''
-examine the connection of drugs to the knockdown of their expected gene target

Larson Hogstrom, 4/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.util.progress as update
from matplotlib import cm

wkdir = '/xchip/cogs/projects/target_id/drug_gene_connections_20Mar2014'

dFile = '/xchip/cogs/projects/target_id/drug_gene_connections_20Mar2014/drug_gene_pairs_02-2014.txt'
targetSheet = pd.read_csv(dFile,sep='\t')

### load summly matrix
### use lass matched matrix
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_n7147x7147.gctx'
matrixType = 'rnkpt_matched_lass'
gt = gct.GCT()
gt.read(mtrxSummly)
summFrm = gt.frame
# load meta-data
summFrm.index = gt.get_column_meta('pert_iname')
summFrm.columns = gt.get_column_meta('pert_id')
pert_type = gt.get_column_meta('pert_type')
pert_iname = gt.get_column_meta('pert_iname')
pert_id = gt.get_column_meta('pert_id')
isSh = [x == 'trt_sh.cgs' for x in pert_type]
isCp = [x == 'trt_cp' for x in pert_type]

### make CGS x cp matrix
cgsFrm = summFrm.ix[isSh,isCp]
cgsRank = cgsFrm.rank(ascending=False,axis=0)

### which of the drug-gene pairs are in summly space
geneSummly = targetSheet.gene.isin(pert_iname)
drugSummly = targetSheet.pert_id.isin(pert_id)
inSummly = geneSummly & drugSummly
summlySheet = targetSheet[inSummly]
outF = wkdir + '/drug_gene_pairs_in_summly_space.txt'
summlySheet.to_csv(outF,index=False,sep='\t')

### loop through drug-gene pairs to find rank
for r in summlySheet.iterrows():
    r[0]





# drug x CGS rank

