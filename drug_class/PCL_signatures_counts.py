'''
-Go through the compounds that make up various PCLs
-Count how many sigatures there are of each compound and the distribution
of cell lines
-identify 'weak' compounds with less signatures 

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

# PCL annotation file
# gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140221/cliques.gmt' #most up-to date drug groups
gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140213/cliques.gmt'
# gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140402/cliques.gmt'

wkdir = '/xchip/cogs/projects/pharm_class/lhwork/pcl_20140213_stats'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

cliqueGMT = gmt.read(gFile)
cliqFrm = pd.DataFrame(cliqueGMT)
# unstack nested list
cliqMemberLong = [item for sublist in cliqFrm.sig.values for item in sublist]
cliqMemb = list(set(cliqMemberLong))

# load summly matrix
summMtrx = '/xchip/cogs/projects/connectivity/summly/matrices/matched_mrp4_n7147x7147.gctx'
gt = gct.GCT()
gt.read(summMtrx)
summFrm = gt.frame

### get info on drug signatures
MC = mu.CMapMongo()
pertInfo = MC.find({'pert_id':{'$in':cliqMemb}},{'sig_id':True,'cell_id':True,'pert_id':True,'pert_iname':True},toDataFrame=True)

# tabulate signature stats
pertGrped = pertInfo.groupby('pert_id')
nDrugs = len(pertGrped.groups)
Zs = np.zeros((nDrugs,5))
sig_counts = pd.DataFrame(Zs,index=pertGrped.groups.keys(),columns=['pert_iname','n_cells','n_signatures','n_PCLs','PCLs'])
sig_counts.index.name = 'pert_id'
for xt in pertGrped:
    brd = xt[0]
    grp = xt[1]
    sig_counts.ix[brd,'pert_iname'] = grp.pert_iname.values[0]
    cell_set = set(grp.cell_id)
    sig_counts.ix[brd,'n_cells'] = len(cell_set)
    sig_counts.ix[brd,'n_signatures'] = grp.shape[0]
    # check which PCLS the compound belongs to
    def is_in(x):
        return brd in x
    inPCL = cliqFrm[cliqFrm.sig.apply(is_in)]
    sig_counts.ix[brd,'n_PCLs'] = inPCL.shape[0]
    sig_counts.ix[brd,'PCLs'] = list(inPCL.desc.values)
sig_counts['is_summly'] = sig_counts.index.isin(summFrm.index)
outF = os.path.join(wkdir, 'PCL_signature_counts.txt')
sig_counts.to_csv(outF, sep='\t')

plt.hist(sig_counts.n_signatures,30,range=(0,200))
plt.ylabel('number of compounds',fontweight='bold')
plt.xlabel('n_signatures ',fontweight='bold')
plt.title('signature counts for PCL members')
outF = os.path.join(wkdir, 'PCL_n_sigatures.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

plt.hist(sig_counts.n_cells,30)
plt.ylabel('number of compounds',fontweight='bold')
plt.xlabel('n_cells ',fontweight='bold')
plt.title('cell line counts for PCL members')
outF = os.path.join(wkdir, 'PCL_n_cell_lines.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

