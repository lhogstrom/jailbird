'''
- load in Rajiv's PCL analysis (clique)
- look for overlap among groups
- look for specificity of PCL gorupings

Larson Hogstrom, 1/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import cmap
import os
import cmap.io.gmt as gmt
import scipy.cluster

wkdir = '/xchip/cogs/projects/pharm_class/lh_clique_6Jan14'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

#load in clique annotations and matrix
cFile = '/xchip/cogs/sig_tools/sig_cliquescore_tool/sample/cp_clique_n69/clique.gmt'
cliqueGMT = gmt.read(cFile)
cliqFrm = pd.DataFrame(cliqueGMT)

cMedianFile = '/xchip/cogs/sig_tools/sig_cliquescore_tool/sample/cp_clique_n69/clique_median_n69x7147.gctx'
cliqueGCT = gct.GCT(cMedianFile)
cliqueGCT.read()
pIDs = cliqueGCT.get_row_meta('pert_id')
cmFrm = cliqueGCT.frame
cmFrm.index = pIDs

# duplicate the PCLxPCL matrix by correlating PCL vectors
corrMtrx = np.corrcoef(cmFrm.values,rowvar=0)
corrFrm = pd.DataFrame(data=corrMtrx,index=cmFrm.columns,columns=cmFrm.columns)
#perform hierarchical clustering
Y = scipy.cluster.hierarchy.linkage(corrMtrx, method='centroid')
Z = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
cOrder = Z['leaves']
iPCL = corrFrm.index[cOrder]
clustered = corrFrm.reindex(index=iPCL,columns=iPCL)
#make heatmap
plt.imshow(clustered.values,interpolation='nearest')
ytcks = list(clustered.index.values)
plt.xticks(np.arange(len(ytcks)), ytcks,rotation=75)
plt.yticks(np.arange(len(ytcks)),ytcks)
plt.colorbar()
outF = os.path.join(wkdir,'pcl_by_pcl_matrix')
plt.savefig(outF, bbox_inches='tight')
# does this replicate with cmap.analytics.cluster? 

# how much overlap is there among groups? 
# Are there specific compound-compound pairs that are driving intra PCL relatedness?
inameDict = {}
for x in cliqFrm.iterrows():
    inameDict[x[1]['id']] = set(x[1]['sig'])

## find compounds which overlap in multiple 
nInm = len(inameDict)
overlapCount = pd.DataFrame(np.zeros([nInm, nInm]),
    index=inameDict.keys(),
    columns=inameDict.keys())
overlapProp = pd.DataFrame(np.zeros([nInm, nInm]),
    index=inameDict.keys(),
    columns=inameDict.keys())
JaccardIndex = pd.DataFrame(np.zeros([nInm, nInm]),
    index=inameDict.keys(),
    columns=inameDict.keys())
for iname1 in inameDict:
    grp1Set = inameDict[iname1]
    for iname2 in inameDict:
        grp2Set = inameDict[iname2]
        interSect = grp1Set.intersection(grp2Set)
        unioN = grp1Set.union(grp2Set)
        nIntersect = len(interSect)
        nUnion = len(unioN)
        overlapCount.ix[iname1,iname2] = nIntersect #overlap counts
        propIntersect = nIntersect/float(len(grp1Set)) #proportion of overlap
        overlapProp.ix[iname1,iname2] = propIntersect
        jaccardI = nIntersect/float(nUnion)
        JaccardIndex.ix[iname1,iname2] = jaccardI
outF = wkdir + '/PCL_overlap_proportion_matrix.txt'
overlapProp.to_csv(outF,sep='\t',index=True,header=True)
outF = wkdir + '/PCL_overlap_count_matrix.txt'
overlapCount.to_csv(outF,sep='\t',index=True,header=True)
outF = wkdir + '/PCL_Jaccard_index_matrix.txt'
JaccardIndex.to_csv(outF,sep='\t',index=True,header=True)
## what are the top overlaps between groups?
upperFrm = overlapProp.copy()
np.fill_diagonal(upperFrm.values, np.nan)
overlapSer = upperFrm.unstack()
overlapSer = overlapSer[~overlapSer.isnull()] #remove nulls 
overlapSer.sort(ascending=False)
overlapSer.name = 'percent_overlap_among_pair'
overlapSer.index.name = 'PCL_pairs'
outF = wkdir +  '/PCL_overlap_proportion_list.txt'
overlapSer.to_csv(outF,sep='\t',index=True,header=True)

### make heatmap of Jaccard Index Matrix
fig = plt.figure(figsize=(20, 20), dpi=50)
plt.imshow(JaccardIndex,
    interpolation='nearest',
    cmap=cm.gray_r)
# ytcks = list(JaccardIndex.index.values)
ytcks = [x[:30] for x in JaccardIndex.index.values]
xtcks = ytcks
# plt.xticks(np.arange(len(xtcks)), xtcks,rotation=75)
plt.xticks(np.arange(len(xtcks)),xtcks,rotation=90)
plt.yticks(np.arange(len(ytcks)),ytcks)
plt.colorbar()
outF = os.path.join(wkdir,'jaccard_index_matrix.png')
# plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()
