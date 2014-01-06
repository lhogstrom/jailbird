'''
- load in Rajiv's PCL analysis (clique)
- look for overlap among groups
- look for specificity of PCL gorupings

Larson Hogstrom, 1/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import cmap
import os
import cmap.io.gmt as gmt

cFile = '/xchip/cogs/sig_tools/sig_cliquescore_tool/sample/cp_clique_n69/clique.gmt'
cliqueGMT = gmt.read(cFile)
cliqFrm = pd.DataFrame(cliqueGMT)

cMedianFile = '/xchip/cogs/sig_tools/sig_cliquescore_tool/sample/cp_clique_n69/clique_median_n69x7147.gctx'
cliqueGCT = gct.GCT(cMedianFile)
cliqueGCT.read()
pIDs = cliqueGCT.get_row_meta('pert_id')
cmFrm = cliqueGCT.frame
cmFrm.index = pIDs



# try to duplicate the PCLxPCL matrix by correlating PCL vectors
corrMtrx = np.corrcoef(cmFrm.values,rowvar=0)
corrFrm = pd.DataFrame(data=corrMtrx,index=cmFrm.columns,columns=cmFrm.columns)


#perform hierarchical clustering
import scipy.cluster
pDist = scipy.spatial.distance.pdist(corrMtrx)
z = scipy.cluster.hierarchy.single(pDist)


Y = scipy.cluster.hierarchy.linkage(corrMtrx, method='centroid')
Z = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
cOrder = Z['leaves']
iPCL = corrFrm.index[cOrder]
clustered = corrFrm.reindex(index=iPCL,columns=iPCL)

#make heatmap
plt.imshow(clustered.values,interpolation='nearest')
ytcks = list(clustered.index.values)
plt.xticks(np.arange(len(yticks)), ytcks,rotation=75)
plt.yticks(np.arange(len(yticks)),ytcks)
plt.colorbar()


plt.title('sum_score')



d = scipy.cluster.hierarchy.distance.pdist(corrMtrx)   # vector of (100 choose 2) pairwise distances
# L = scipy.cluster.hierarchy.linkage(d, method='complete')
L = scipy.cluster.hierarchy.linkage(d, method='complete')
ind = scipy.cluster.hierarchy.fcluster(L, 0.5*d.max(), 'distance')



#compute the histogram of the linkage distances and find the first empty bin. Use that
#bin to set the cutoff for cluster separation
zdist = [z[i,2] for i in range(len(z))]
n,bins,patches = plt.hist(zdist,bins=bins) #@UnusedVariable
plt.clf()
try:
    distance_cutoff = bins[np.where(n == 0)][0]
except IndexError:
    distance_cutoff = np.max(bins)

#cut the linkage hierarchy at the computed distance cutoff and assign cluster membership
t = scipy.cluster.hierarchy.fcluster(z, distance_cutoff, criterion='distance')
#build a 2D list in which each entry is a list of data points that fall in each cluster
num_clusters = np.max(t)
clusters = []
cluster_data = []
for i in range(num_clusters):
    clusters.append([labels[j] for j in range(len(t)) if t[j] == i+1])
    cluster_data.append([input_sequence[j] for j in range(len(t)) if t[j] == i+1])

# #return the list of clusters
# return (clusters,cluster_data)
