#!/bin/py 

import cmap.io.gct as gct
from pylab import plot,show
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq

# data generation
data = vstack((rand(150,2) + array([.5,.5]),rand(150,2)))

# computing K-Means with K = 2 (2 clusters)
centroids,_ = kmeans(data,2)
# assign each sample to a cluster
idx,_ = vq(data,centroids)

# some plotting using numpy's logical indexing
plot(data[idx==0,0],data[idx==0,1],'ob',
     data[idx==1,0],data[idx==1,1],'or')
plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
show()

### try clustering some expression data
fname = '/xchip/obelix/pod/brew/pc/NMH002_NEU_6H/by_pert_id/NMH002_NEU_6H_COMPZ.MODZ_SCORE_LM_n336x978.gctx'
db = gct.GCT()
db.read(fname)

x = db.matrix
#to transpose or not to transpose - that is the question
xt = x.T #The columns are the features seen during each observation.
#probably don't need to use vq 'whiten' because the data is in z-score space
centroids,_ = kmeans(xt,13)
idx,_ = vq(xt,centroids) #these are the cluster asignments for each

iByCluster = numpy.argsort(idx)
clustOrder = idx[iByCluster]
xClustered = xt[iByCluster,:]

#work on showing cluster labels
fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
ax = axs[0]
ax.imshow(xt,cmap=cmap)
ax.set_title('before clustering')
ax = axs[1]
ax.imshow(xClustered,cmap=cmap)
ax.set_title('after k means clustering')
# fig.suptitle('')
plt.show()



# ind2 = ind2[:,idx2]
# plt.imshow(xt,cmap=cmap)
# plt.imshow(xClustered,cmap=cmap)
