'''
-examine NMF results across cell lines
-build benchmarks

Larson Hogstrom, 12/2013
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import os
import cmap.analytics.summly_null as SN
from statsmodels.distributions import ECDF
import cmap.io.gmt as gmt
import cmap.util.progress as update
from matplotlib import cm

wkdir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/NMF_benchmark_development'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

### Analyze W matrix ###
# wkdir = '/xchip/cogs/projects/NMF/clique_n69_all_cell_lines'
Hfile = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c9_LM/clique_compound_classes_n585x978.H.k9.gct'
WFile = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c9_LM/clique_compound_classes_n585x978.W.k9.gct'
aFile = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c9_LM/clique_compound_classes.v3.txt'
# Hmtrx = pd.io.parsers.read_csv(Hfile,sep='\t',index_col=0) #,header=True
Hmtrx = pd.read_csv(Hfile,sep='\t',skiprows=[0,1],index_col=0) #,header=True
Hmtrx = Hmtrx.drop('Description',1)
Hmtrx = Hmtrx.T
# Wmtrx = pd.io.parsers.read_csv(WFile,sep='\t',index_col=0) #,header=True
anntFrm = pd.read_csv(aFile,sep='\t',header=False,index_col=0)
anntFrm.columns = reducedSigFrm.columns

### load in clique annotations and matrix
cFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140221/cliques.gmt'
cliqueGMT = gmt.read(cFile)
cliqFrm = pd.DataFrame(cliqueGMT)

# groupSer = pd.Series(index=anntFrm.index,data=anntFrm.pcl_name)
# if (groupSer.index == Hmtrx.index).all():
#     iZip = zip(*[groupSer.values,groupSer.index.values])
#     mCol = pd.MultiIndex.from_tuples(iZip, names=['pcl_name','sig_id'])
#     Hmtrx.index = mCol
# pclGrped = Hmtrx.groupby(level='pcl_name')
graphDir = wkdir + '/graphs_uniform_max_sort'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
# flH = np.float64(Hmtrx.values)
maxVal = Hmtrx.max(axis=1).max()
for r in cliqFrm.iterrows():
    grp = r[1]['id']
    brds = r[1]['sig']
    anntMtch = anntFrm[anntFrm.pert_id.isin(brds)]
    grpH = Hmtrx.reindex(anntMtch.index)
    meanVec = grpH.describe().ix['mean']
    ### correlate each observed sigature with mean vector
    # corrMtrx = np.corrcoef(meanVec,grpH)
    # corrVec = corrMtrx[:-1,-1]
    # iSort = np.argsort(corrVec)
    # grpH = grpH.ix[iSort,:] # sort acording to corr with mean
    ### take top three components - order acording to their strenght
    iTop3 = meanVec.order(ascending=False).index[:3]
    sortedTop = grpH.ix[:,iTop3].sort()
    topSum = sortedTop.sum(axis=1).order(ascending=False)
    grpH = grpH.ix[topSum.index,:] # sort acording to corr with mean
    Hfloat = np.float64(grpH.values)
    fig = plt.figure(figsize=(40, 20), dpi=50)
    plt.imshow(Hfloat,
        interpolation='nearest',
        cmap=cm.gray_r,
        vmax=maxVal)
    ytcks = list(grpH.index)
    xtcks = list(grpH.columns)
    plt.xticks(np.arange(len(xtcks)), xtcks,rotation=75)
    plt.yticks(np.arange(len(ytcks)),ytcks)
    plt.colorbar()
    grpMod = grpMod = ''.join(e for e in grp if e.isalnum())
    outF = os.path.join(graphDir,grpMod+'.png')
    plt.savefig(outF, bbox_inches='tight')
    plt.close()

### load NMF projection results
gFile = '/xchip/cogs/projects/NMF/MCF7_comp_annot_to_CCLE_space2/c_annotc1/c1_vs_ACHILLES_Comp_annot.v1.pdf.gct'
gt = gct.GCT()
gt.read(gFile)

