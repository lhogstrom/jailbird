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
isSh = [x == 'trt_sh.cgs' for x in pert_type] # trt_oe
isCp = [x == 'trt_cp' for x in pert_type]
pert_iname_dict = {}
for ix, x in enumerate(pert_id):
    pert_iname_dict[x] = pert_iname[ix]

### make CGS x cp matrix
cgsFrm = summFrm.ix[isSh,isCp]
cgsRank = cgsFrm.rank(ascending=False,axis=0)

### which of the drug-gene pairs are in summly space
geneSummly = targetSheet.gene.isin(cgsRank.index)
drugSummly = targetSheet.pert_id.isin(cgsRank.columns)
inSummly = geneSummly & drugSummly
summlySheet = targetSheet[inSummly]
outF = wkdir + '/drug_gene_pairs_in_summly_space.txt'
summlySheet.to_csv(outF,index=False,sep='\t')

### loop through drug-gene pairs to find rank
zer = np.zeros([len(summlySheet),5])
dgIndex = summlySheet.pert_id + ':' + summlySheet.gene
resFrame = pd.DataFrame(zer,index=dgIndex,columns=['drug_id','drug_iname','gene','connection_rank','mean_rank_pt_4'])
for r in summlySheet.iterrows():
    gene = r[1][0]
    pID = r[1][1]
    di = pID + ':' + gene
    pIname = pert_iname_dict[pID]
    resFrame.ix[di,'drug_id'] = pID
    resFrame.ix[di,'drug_iname'] = pIname
    resFrame.ix[di,'gene'] = gene
    resFrame.ix[di,'connection_rank'] = cgsRank.ix[gene,pID]
    resFrame.ix[di,'mean_rank_pt_4'] = cgsFrm.ix[gene,pID]
resFrame = resFrame.sort('connection_rank')
outF = wkdir + '/expected_drug_gene_connection_ranks.txt'
resFrame.to_csv(outF,index=False,sep='\t')

### plot ranks
fig = plt.figure(1, figsize=(8,8))
plt.plot(resFrame.connection_rank,'.')
outF = wkdir + '/expected_drug_gene_connection_ranks.png'
plt.ylabel('rank out of 3797 genes')
plt.xlabel('summly drug-gene connections')
plt.ylim((1,cgsRank.shape[0]))
plt.xlim((1,resFrame.shape[0]))
plt.title('expected drug connection to gene KD')
plt.savefig(outF, bbox_inches='tight')
plt.close()
