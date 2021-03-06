'''
-examine the is_gold status of compounds across cell lines

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

wkdir = '/xchip/obelix/pod/non_public/ROC/lhwork/is_gold_evaluation'

# load summary file
sFile = '/xchip/obelix/pod/non_public/ROC/build/with_dmso/siginfo.txt'
si = pd.read_csv(sFile,sep='\t')

### make compound x cell line matrix - percent gold
uCp = set(si.pert_id)
uCell = set(si.cell_id)
nCp = len(uCp)
nCell = len(uCell)
mZero = np.zeros([nCp,nCell])
percent_gold = pd.DataFrame(mZero, index=list(uCp),columns=list(uCell))
gold_counts = pd.DataFrame(mZero, index=list(uCp),columns=list(uCell))
pertGrped = si.groupby(['pert_id','cell_id'])
si.is_gold = si.is_gold.replace(-666,np.nan)
for grp in pertGrped.groups:
    brd = grp[0]
    cell = grp[1]
    igroup = pertGrped.groups[grp]
    grpFrm = si.ix[igroup,:]`
    gold_percent = sum(grpFrm['is_gold'])/float(grpFrm.shape[0])
    percent_gold.ix[brd,cell] = gold_percent
    gold_counts.ix[brd,cell] = sum(grpFrm['is_gold'])
outFile = wkdir + '/compound_by_cell_line_pecent_is_gold.txt'
percent_gold.to_csv(outFile, sep='\t')
# 1. For each pert_id calculate number of cell lines that have >= 50% signatures that are gold.
mostly_gold = percent_gold >= .5
mostly_gold_gold.to_csv(outFile, sep='\t')
goldSum = mostly_gold.sum(axis=1)
goldSum.name = 'n_cell_lines_majority_gold'
outFile = wkdir + '/compound_cell_line_counts_with_majority_gold_signatures.txt'


# 2. Compute gold signatures by cell type
# gold counts per cell line
gold_counts_per_cell = gold_counts.sum(axis=1)
gold_counts_per_cell.name = 'gold_counts_per_cell_line'
outFile = wkdir + '/gold_counts_per_cell_line.txt'
gold_counts_per_cell.to_csv(outFile,sep='\t')

# 3. Examine well DMSO stats
dmsoFrm = si[si.pert_iname == 'DMSO']
dmsoFrm.is_gold = dmsoFrm.is_gold.replace(-666,np.nan)
dmso_gold_count = np.nansum(dmsoFrm.is_gold)

dSS = si[si.pert_iname == 'DMSO'].distil_ss
dR1 = si[si.pert_type == 'trt_cp'].distil_ss
h2 = plt.hist(dSS,30,color='r',range=[0,20],label='DMSO',alpha=.5,normed=True)
h1 = plt.hist(dR1.values,30,color='b',range=[0,20],label=['Roche'],alpha=.6,normed=True)
plt.legend()
plt.xlabel('distill_ss')
plt.title('signature strength')
outFile = wkdir + '/ss_hist.png'
plt.savefig(outFile)
plt.close()

dSS = si[si.pert_iname == 'DMSO'].distil_cc_q75
dR1 = si[si.pert_type == 'trt_cp'].distil_cc_q75
h2 = plt.hist(dSS,30,color='r',range=[-.2,1],label='DMSO',alpha=.5,normed=True)
h1 = plt.hist(dR1.values,30,color='b',range=[-.2,1],label=['Roche'],alpha=.6,normed=True)
plt.legend()
plt.xlabel('distil_cc_q75')
plt.title('replicate correlation')
outFile = wkdir + '/cc_hist.png'
plt.savefig(outFile)
plt.close()


# hist of dmso vs is_gold
dmsoFrm.distil_ss
dmsoFrm.distil_cc_q75




# Use doese info - Are gold doses close toghether? 