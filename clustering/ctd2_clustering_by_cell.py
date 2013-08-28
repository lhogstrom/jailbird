#! /usr/bin/env python
'''
take every ctd2 signature from one cell line
put them in a big matrix and cluster them
repeat for each cell line

which compounds stay together most closely and consistantly? - Average distance of drug-drug pairs
does this reflect class labels?

histogram:
cluster DMSOs is the average pairwise correlation higher or lower?
'''
import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd
import matplotlib
from cmap.analytics.pert_explorer import PertExplorer
from cmap.analytics.pert_explorer import BaseLine

from cmap.analytics.cluster import HClust
import cmap
from os import path
from matplotlib import cm

wkdir = '/xchip/cogs/projects/target_id/ctd2_cp_clustering_28Aug2'
if not os.path.exists(wkdir):
    os.makedirs(wkdir)

#what are the cps that were tested in all 50 cell lines?
CM = mu.CMapMongo()    
qr = CM.find({'pert_type':'trt_cp','sig_id' : {'$regex' : 'CPC006'}},{'cell_id':True,'pert_id':True},toDataFrame = True)
ctdCells = list(set(qr['cell_id']))
ctdCells.sort()
ctdCps = list(set(qr['pert_id']))
ctdCps.sort()

cell = 'LOVO'
cpQ = CM.find({'pert_type':'trt_cp','cell_id':cell,'sig_id' : {'$regex' : 'CPC006'}},{'sig_id':True,'pert_id':True},toDataFrame = True)
lim1 = CM.find({'pert_id':{'$in':list(ctdCps)},'cell_id':cell,'sig_id' : {'$regex' : 'CPC006'}},{'sig_id':True,'pert_id':True},limit=1,toDataFrame = True)


### are there any repeated perts? 
#cps from non-core cell line - LOVO
# cpQ = CM.find({'pert_type':'trt_cp','cell_id':'SKM1','sig_id' : {'$regex' : 'CPC006'}},{'sig_id':True,'pert_id':True},toDataFrame = True)
# ctdCps = list(cpQ['pert_id'])

# grpID = cpQ.groupby('pert_id')
# for sig in grpID['sig_id']:
#     print len(sig)

# cnts = grpID.count()
# repeats = cnts[cnts['sig_id'] >1]
# rpt = repeats.index.values
# for x in rpt:
#     print cpQ['sig_id'][cpQ['pert_id']==x]
# cpQ['sig_id'][rpt]

# DMSOs in each cell line:
dmsoCounts = []
for cell in ctdCells:
    dmso1 = CM.find({'pert_id':'DMSO','cell_id':cell,'sig_id' : {'$regex' : 'CPC006'}},{'sig_id':True})
    n_dmso = len(dmso1)
    dmsoCounts.append(n_dmso)


### get one sig_id for each compound and cell line
SigDict = {}
skippedDict = {}
for cell in ctdCells:
    cellSigs = []
    cellSkipped = []
    for cp in ctdCps:
        lim1 = CM.find({'pert_id':cp,'cell_id':cell,'sig_id' : {'$regex' : 'CPC006'}},{'sig_id':True},limit=1)
        if lim1:
            cellSigs.append(lim1[0])
        else:
            cellSkipped.append(cp)
    SigDict[cell] = cellSigs
    skippedDict[cell] = cellSkipped


### perform clustering
metric = 'spearman'
n_ctd = len(ctdCps)
cellDistMtrx= np.zeros((len(ctdCells),n_ctd,n_ctd))
scoreDict = {}
for icell,cell in enumerate(SigDict.keys()):
    out = wkdir + '/' + cell
    if not os.path.exists(out):
        os.mkdir(out)
    sigs = SigDict[cell]
    query = {'sig_id':{'$in':sigs}}
    PE = PertExplorer(pert_id = None,
                        metric = metric,
                        query = query,
                        out = out)
    ### cluster the data, store it on the PertExplorer instance
    PE.hclust = HClust(PE.score.copy(), PE.annots.copy(), out = PE.out,
                         symmfun = np.max, cutoff = 0.5)
    PE.hclust.cluster()
    PE.hclust.get_cluster_order()
    clust_order = PE.hclust.cluster_order.order().index
    # PE.hclust.draw_dendrogram(orientation = 'right',showfig = False)
    # PE.hclust.save_dendrogram()
    ### add cluster assignment to annotations
    PE.annots['hclust_assignment'] = PE.hclust.cluster_assignment
    ### draw heatmap, ordered by the clustering; save
    PE.draw_heatmap(order = clust_order,
                      params = {'ticklabel_fontsize': 8},
                      vmin = -1, vmax = 1,
                      cmap = cm.RdBu_r,
                      col_label = None,
                      row_label = 'pert_iname',
                      annots_top = [('hclust_assignment', False)],
                      annots_right = ['pert_dose',
                                      ('distil_cc_q75', {'vmin' : 0, 'vmax' : 1}),
                                      ('distil_ss', {'vmin' : 0, 'vmax' : 10})],
                      show_title = False,
                      title = cell + ' Clustering ctd2', 
                      showfig = False) #how to turn off row_labels
                      # annots_bottom = [('cell_lineage', True)],
    PE.heatmap.save(format = 'png')
    plt.close()
    ### draw baseline data; save
    # PE.showBaselines(order = PE.annots.cell_id[clust_order].values,
    #                    clust_assignment = PE.hclust.cluster_assignment.order().values)
    # PE.saveBaselines(format = 'png')
    ### write html file and annotations
    # PE.write_html(extra_figs = [(PE.hclust.dendro_fname, 'dendrogram')])
    PE.write_annotations(encoding = 'latin-1') 
    ### calculate distance consistancy
    idClust = PE.heatmap._get_labels('pert_id')
    clustSer = pd.Series(idClust)
    if len(idClust) != len(set(idClust)):
        print "warning: every pert must be unique"
    distanceMtrx = np.zeros((n_ctd,n_ctd))
    distanceDict = {}
    #place cluster distance measures in a matrix
    for ipert1,pert1 in enumerate(ctdCps):
        if pert1 not in PE.annots['pert_id'].values:
            distanceMtrx[ipert1,:] = np.nan
            distanceMtrx[:,ipert1] = np.nan
            continue
        iloc1 = clustSer[clustSer == pert1].index
        loc1 = iloc1[0]
        for ipert2,pert2 in enumerate(ctdCps):
            if pert2 not in PE.annots['pert_id'].values:
                continue
            if pert1 == pert2:
                continue
            else:
                iloc2 = clustSer[clustSer == pert2].index
                loc2 = iloc2[0]
                distance = abs(loc1-loc2)
                distanceMtrx[ipert1,ipert2] = distance
    cellDistMtrx[icell,:,:] = distanceMtrx
    scoreDict[cell] = PE

### calculate average correlations across cell lines
cellScoreMtrx= np.zeros((len(ctdCells),n_ctd,n_ctd))
for icell,cell in enumerate(scoreDict.keys()):
    #reindex to brd
    score = scoreDict[cell].score
    sigs = scoreDict[cell].annots.index
    # do index match annotations
    # score.index == sigs
    # score.columns == sigs
    pIDs = scoreDict[cell].annots['pert_id']
    idScore = score.copy()# replace the index with pert_ids (from sig_ids)
    idScore.index = pIDs
    idScore.columns = pIDs
    idScore = idScore.reindex(index=ctdCps,columns=ctdCps) #reindex to same size for every cell line
    cellScoreMtrx[icell,:,:] = idScore.as_matrix()
### calculate average scores across cell line
csm = cellScoreMtrx[:5,:,:]
mskCsm = np.ma.masked_array(csm,np.isnan(csm))
scoreMeans = np.mean(mskCsm,axis=0)
triScores = np.triu(scoreMeans.data)





### calculate average clustering distances
cdm = cellDistMtrx[:5,:,:]
np.mean(cdm,axis=0)
mskCdm = np.ma.masked_array(cdm,np.isnan(cdm))
distMeans = np.mean(mskCdm,axis=0)
triMeans = np.triu(distMeans.data)
# make dict/Series of all pairwise distances
avDistDict = {}
for ipert1,pert1 in enumerate(ctdCps):
    for ipert2,pert2 in enumerate(ctdCps):
        if triMeans[ipert1,ipert2] == 0:
            continue
        else:
            pairStr = ctdCps[ipert1] + ':' + ctdCps[ipert2]
            avDistDict[pairStr] = triMeans[ipert1,ipert2]
avDistSer = pd.Series(avDistDict)
avDistSer.sort()
# plot the average pairwise distances
plt.plot(avDistSer)
plt.show()

### make pd frame of pairwise names, dist averages, corr averages


# ROC

#to do:
#average correlation across cell lines
#clustering of DMSOs 



### scratch
il = np.triu_indices(len(distMeans),k=0)
upMeans = distMeans.data.copy()
upMeans[il] = np.nan
randFlat = upMeans.flatten()
uniqRand = upMeans[~np.isnan(upMeans)]

# calculate average distance for each drug-drug connection
for ipert1,pert1 in enumerate(ctdCps):
    for ipert2,pert2 in enumerate(ctdCps):
        distArray = cdm[:,ipert1,ipert2]
        mDistArray = np.ma.masked_array(distArray,np.isnan(distArray))
        avDist = np.mean(mDistArray)    


        mdat = np.ma.masked_array(dat,np.isnan(dat))
mm = np.mean(mdat,axis=1)
