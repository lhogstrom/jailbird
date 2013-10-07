#!/bin/py

### modified from - hierarchical_clustering.py
#2005-2012 J. David Gladstone Institutes, San Francisco California

import matplotlib.pyplot as pylab
import matplotlib.pyplot as plt
from matplotlib import mpl
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy
import string
import time
import sys, os
import getopt

row_method = 'average'
column_method = 'single'
row_metric = 'cityblock' #cosine
column_metric = 'euclidean'
color_gradient = 'red_white_blue'
cmap=pylab.cm.bwr

filename = '/xchip/cogs/projects/NMH/cfwork/hierarch_clustering/NMH001_NEU_6H_COMPZ_tabTxt.txt'

def importData(filename):
    start_time = time.time()
    matrix=[]
    row_header=[]
    first_row=True

    if '/' in filename:
        dataset_name = string.split(filename,'/')[-1][:-4]
    else:
        dataset_name = string.split(filename,'\\')[-1][:-4]
        
    for line in open(filename,'rU').xreadlines():         
        t = string.split(line[:-1],'\t') ### remove end-of-line character - file is tab-delimited
        if first_row:
            column_header = t[1:]
            first_row=False
        else:
            if ' ' not in t and '' not in t: ### Occurs for rows with missing data
                s = map(float,t[1:])
                if (abs(max(s)-min(s)))>0:
                    matrix.append(s)
                    row_header.append(t[0])
            
    time_diff = str(round(time.time()-start_time,1))
    try:
        print '\n%d rows and %d columns imported for %s in %s seconds...' % (len(matrix),len(column_header),dataset_name,time_diff)
    except Exception:
        print 'No data in input file.'; force_error
    return numpy.array(matrix), column_header, row_header

#import data
matrix, column_header, row_header = importData(filename)
x = matrix


### Scale the Matplotlib window size
default_window_hight = 8.5
default_window_width = 12
fig = pylab.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
color_bar_w = 0.015 ### Sufficient size to show

### Scale the max and min colors so that 0 is white/black
vmin=x.min()
vmax=x.max()
vmax = max([vmax,abs(vmin)])
vmin = vmax*-1
norm = mpl.colors.Normalize(vmin/2, vmax/2) ### adjust the max and min to scale these colors

#heatmap code
# Compute and plot top dendrogram

start_time = time.time()
d2 = dist.pdist(x.T) #Computes the pairwise distances between m original observations in n-dimensional space.
D2 = dist.squareform(d2) #Converts a vector-form distance vector to a square-form distance matrix
Y2 = sch.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
### linkage output - array - The hierarchical clustering encoded as a linkage matrix.
Z2 = sch.dendrogram(Y2) #
ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
#the treshold is set to t=0.7*max(Y2[:,2]) - why this value?
time_diff = str(round(time.time()-start_time,1))
print 'Column clustering completed in %s seconds' % time_diff

start_time = time.time()
d1 = dist.pdist(x) #Computes the pairwise distances between m original observations in n-dimensional space.
D1 = dist.squareform(d1) #Converts a vector-form distance vector to a square-form distance matrix
Y1 = sch.linkage(D1, method=row_method, metric=row_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
### linkage output - array - The hierarchical clustering encoded as a linkage matrix.
Z1 = sch.dendrogram(Y1) #
ind2 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
#the treshold is set to t=0.7*max(Y2[:,2]) - why this value?
time_diff = str(round(time.time()-start_time,1))
print 'Column clustering completed in %s seconds' % time_diff

###apply the clustering index to the matrix of data
xt = x
idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
xt = xt[:,idx2]
ind2 = ind2[:,idx2]
plt.imshow(xt,cmap=cmap)


