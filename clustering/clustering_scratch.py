'''
scratch code for clustering analysis
'''

import cmap.analytics.cluster as clst
m1 = np.random.rand(6,3)
a1 = np.random.rand(6,2)
l1 = ['a','b','c','d','e','f']
m1 = pd.DataFrame(m1)
m1.index = l1
# m1.columns = l1
annots1 = pd.DataFrame(l1,index=l1)
cO = clst.APClust(m1,annots1,convert_to_distance=True)
cl = clst.CMapClust()
sl = clst.single_linkage_hist_cluster(m1)

#perform hierarchical clustering on 
import scipy.cluster
Y = scipy.cluster.hierarchy.linkage(fpFrame, method='centroid')
Z = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
cOrder = Z['leaves']
iPCL = fpFrame.index[cOrder]
clustered = fpFrm.reindex(index=iPCL,columns=iPCL)
test 2
m1 = np.random.rand(10,10)
m1 = pd.DataFrame(m1)
Y = scipy.cluster.hierarchy.linkage(m1, method='centroid')    
Z = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
cOrder = Z['leaves']
iPCL = m1.index[cOrder]
clustered = m1.reindex(index=iPCL,columns=iPCL)    
#use cmap clustering tool
