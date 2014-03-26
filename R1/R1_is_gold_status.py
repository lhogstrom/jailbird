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
pertGrped = si.groupby(['pert_id','cell_id'])
for grp in pertGrped.groups:
    brd = grp[0]
    cell = grp[1]
    igroup = pertGrped.groups[grp]
    grpFrm = si.ix[igroup,:]
    gold_percent = sum(grpFrm['is_gold'])/float(grpFrm.shape[0])
    percent_gold.ix[brd,cell] = gold_percent
outFile = wkdir + 'compound_by_cell_line_pecent_is_gold.txt'
# 1. For each pert_id calculate number of cell lines that have >= 50% signatures that are gold.
mostly_gold = percent_gold > .5

# 2. Compute gold signatures by cell type
# gold counts per cell line


# 3. Examine well DMSO stats
dmsoFrm = si[si.pert_iname == 'DMSO']
dmsoFrm.is_gold = dmsoFrm.is_gold.replace(-666,np.nan)
dmso_gold_count = np.nansum(dmsoFrm.is_gold)

# hist of dmso vs is_gold
dmsoFrm.distil_ss
dmsoFrm.distil_cc_q75




# Use doese info - Are gold doses close toghether? 