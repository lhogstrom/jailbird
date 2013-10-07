#!/bin/py

### modified from - hierarchical_clustering.py
#2005-2012 J. David Gladstone Institutes, San Francisco California

import matplotlib.pyplot as pylab
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
    
## calculate positions for all elements
# ax1, placement of dendrogram 1, on the left of the heatmap
#if row_method != None: w1 = 
[ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.22,0.2,0.6]   ### The second value controls the position of the matrix relative to the bottom of the view
width_between_ax1_axr = 0.004
height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix

# axr, placement of row side colorbar
[axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
axr_x = ax1_x + ax1_w + width_between_ax1_axr
axr_y = ax1_y; axr_h = ax1_h
width_between_axr_axm = 0.004

# axc, placement of column side colorbar
[axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing
axc_x = axr_x + axr_w + width_between_axr_axm
axc_y = ax1_y + ax1_h + height_between_ax1_axc
height_between_axc_ax2 = 0.004

# axm, placement of heatmap for the data matrix
[axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
axm_x = axr_x + axr_w + width_between_axr_axm
axm_y = ax1_y; axm_h = ax1_h
axm_w = axc_w

# ax2, placement of dendrogram 2, on the top of the heatmap
[ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.15] ### last one controls hight of the dendrogram
ax2_x = axr_x + axr_w + width_between_axr_axm
ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
ax2_w = axc_w

# axcb - placement of the color legend
[axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.88,0.18,0.09]

#heatmap code
# Compute and plot top dendrogram
if column_method != None:
    start_time = time.time()
    d2 = dist.pdist(x.T) #Computes the pairwise distances between m original observations in n-dimensional space.
    D2 = dist.squareform(d2) #Converts a vector-form distance vector to a square-form distance matrix
    ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)
    Y2 = sch.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
    ### linkage output - array - The hierarchical clustering encoded as a linkage matrix.
    Z2 = sch.dendrogram(Y2) #this plots the dendrogram
    ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
    #the treshold is set to t=0.7*max(Y2[:,2]) - why this value?
    ax2.set_xticks([]) ### Hides ticks
    ax2.set_yticks([])
    time_diff = str(round(time.time()-start_time,1))
    print 'Column clustering completed in %s seconds' % time_diff
else:
    ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data
        
# Compute and plot left dendrogram.
if row_method != None:
    start_time = time.time()
    d1 = dist.pdist(x)
    D1 = dist.squareform(d1)  # full matrix
    ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True) # frame_on may be False
    Y1 = sch.linkage(D1, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
    Z1 = sch.dendrogram(Y1, orientation='right')
    ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
    ax1.set_xticks([]) ### Hides ticks
    ax1.set_yticks([])
    time_diff = str(round(time.time()-start_time,1))
    print 'Row clustering completed in %s seconds' % time_diff
else:
    ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data
    
# Plot distance matrix.
axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
xt = x
if column_method != None:
    idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
    xt = xt[:,idx2]
    ind2 = ind2[:,idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
if row_method != None:
    idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
    xt = xt[idx1,:]   # xt is transformed x
    ind1 = ind1[idx1,:] ### reorder the flat cluster to match the order of the leaves the dendrogram
### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black
axm.set_xticks([]) ### Hides x-ticks
axm.set_yticks([])

# Add text
new_row_header=[]
new_column_header=[]
for i in range(x.shape[0]):
    if row_method != None:
        if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
            axm.text(x.shape[1]-0.5, i, '  '+row_header[idx1[i]])
        new_row_header.append(row_header[idx1[i]])
    else:
        if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
            axm.text(x.shape[1]-0.5, i, '  '+row_header[i]) ### When not clustering rows
        new_row_header.append(row_header[i])
for i in range(x.shape[1]):
    if column_method != None:
        axm.text(i, -0.9, ' '+column_header[idx2[i]], rotation=270, verticalalignment="top") # rotation could also be degrees
        new_column_header.append(column_header[idx2[i]])
    else: ### When not clustering columns
        axm.text(i, -0.9, ' '+column_header[i], rotation=270, verticalalignment="top")
        new_column_header.append(column_header[i])

# Plot colside colors
# axc --> axes for column side colorbar
if column_method != None:
    axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
    cmap_c = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
    dc = numpy.array(ind2, dtype=int)
    dc.shape = (1,len(ind2)) 
    im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
    axc.set_xticks([]) ### Hides ticks
    axc.set_yticks([])

# Plot rowside colors
# axr --> axes for row side colorbar
if row_method != None:
    axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for column side colorbar
    dr = numpy.array(ind1, dtype=int)
    dr.shape = (len(ind1),1)
    #print ind1, len(ind1)
    cmap_r = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
    im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
    axr.set_xticks([]) ### Hides ticks
    axr.set_yticks([])

# Plot color legend
axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
axcb.set_title("colorkey")

if '/' in filename:
    dataset_name = string.split(filename,'/')[-1][:-4]
    root_dir = string.join(string.split(filename,'/')[:-1],'/')+'/'
else:
    dataset_name = string.split(filename,'\\')[-1][:-4]
    root_dir = string.join(string.split(filename,'\\')[:-1],'\\')+'\\'
filename = root_dir+'Clustering-%s-hierarchical_%s_%s.pdf' % (dataset_name,column_metric,row_metric)
cb.set_label("Differential Expression (log2 fold)")
# exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2)

### Render the graphic
if len(row_header)>50 or len(column_header)>50:
    pylab.rcParams['font.size'] = 5
else:
    pylab.rcParams['font.size'] = 8

pylab.savefig(filename)
print 'Exporting:',filename
filename = filename[:-3]+'png'
pylab.savefig(filename, dpi=100) #,dpi=200
pylab.show()
