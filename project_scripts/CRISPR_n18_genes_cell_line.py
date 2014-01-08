'''
-First CRISPR pilot will be 18 genes tested in two cell lines.
-First cell line is A375. Need to pick the second cell line.

'''
import cmap.analytics.pcla as pcla
import numpy as np
import os
import cmap.io.gmt as gmt
import cmap
import pandas as pd
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu

wkdir = '/xchip/cogs/projects/CRISPR/pilot'

cFile = '/xchip/cogs/projects/CRISPR/pilot/CRISPR.txt'
crFrm = pd.read_csv(cFile, sep='\t')
#limit to plates FXA06 & FXA07
# plateKeep = ['FXA006','FXA007'] #also 8-12
plateKeep = ['FXA012','FXA007', 'FXA006', 'FXA010','FXA009','FXA008','FXA011']
plateFrm = crFrm[crFrm.Plate.isin(plateKeep)]
gSet = set(plateFrm['Target Gene'])


#baseline expression in the core cell lines
#dynamic expression range
mc = mu.MongoContainer()
ci = mc.gene_info.find({'pr_gene_symbol':{'$in':list(gSet)}},
            {'is_expressed':True,'pr_gene_symbol':True},
            toDataFrame=True)
grpedCi = ci.groupby('pr_gene_symbol')
exprFrm = grpedCi.first()
ciSet = set(ci['pr_gene_symbol'])

#construct expression data in a dataframe
beFrm = pd.DataFrame()
for ix1 in exprFrm.index:
    gene = ix1
    bEx = exprFrm.ix[ix1,'is_expressed']
    bSer = pd.Series(bEx)
    bSer.name = gene
    bFrm = pd.DataFrame(bSer)
    beFrm = pd.concat([beFrm,bFrm],axis=1)
# coreList = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
coreList = ['A375','A549', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # take out HA1E since there is no info for that
# reindex to just core cell lines
coreFrm = beFrm.reindex(coreList)
outF = wkdir + '/crispr_pilot_baseline_exprn_FXA006-12.txt'
coreFrm.to_csv(outF,sep='\t')
plate_expression = coreFrm.sum(axis=1)
