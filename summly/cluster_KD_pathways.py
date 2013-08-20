#! /usr/bin/env python
'GPPA (genomic perturbation pathay analyzer) - cluster genomic perturbation data acording to pathways'

import numpy as np
import os
import cmap.util.mongo_utils as mu
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
from cmap.analytics.pert_explorer import PertExplorer
from cmap.analytics.cluster import HClust
import cmap.io.gmt as gmt
from matplotlib import cm
from cmap.analytics.queryer import Queryer

wkdir = '/xchip/cogs/projects/target_id/KD_pathway_clustering/pert_explorer_mtor_pathway2'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

#pathway annotations from reactome 
pathGMT = '/xchip/cogs/projects/target_id/KD_pathway_clustering/ReactomePathways.gmt'
gmtDict = gmt.read(pathGMT)
pathwayDict = {}
for dict1 in gmtDict:
    pathwayDict[dict1['id']] = dict1['sig']


# 'Cholesterol biosynthesis' # vs statin signature?
# 'p53-Dependent G1 DNA Damage Response'
# 'p53-Dependent G1/S DNA damage checkpoint'
# 'Antigen processing: Ubiquitination & Proteasome degradation'
# 'Regulation of activated PAK-2p34 by proteasome mediated degradation'
# 'mTOR signalling'
# Lessons from the cancer genome pathways:
# 'Signaling by TGF-beta Receptor Complex'
pathName = 'mTOR signalling'
mtorGenes = pathwayDict[pathName]
geneAll = set(mtorGenes)
CM = mu.CMapMongo()
mtorsigs = CM.find({'pert_type':'trt_sh.cgs','pert_iname':{'$in':list(geneAll)}},{'sig_id':True,'pert_iname':True},toDataFrame = True)
mtorcells = CM.find({'pert_type':'trt_sh.cgs','pert_iname':{'$in':list(geneAll)}},{'cell_id':True},toDataFrame = True)
# cellSet = set(mtorcells['cell_id'])
# cellSet.remove()
cell_counts = mtorcells.groupby('cell_id').size()
filt = cell_counts[cell_counts > 2] #filter cell lines with less than 3 pathway signatures in them
cellSet  = set(filt.index)

#load in onco gene list
# inFile = '/xchip/cogs/projects/oncoDome/OncoDome_genes.txt'
# cgsList = []
# with open(inFile,'rt') as f:
#     for string in f:
#         splt = string[:-1]
#         cgsList.append(splt)
# geneAll = set(cgsList)
# # check to see if gene has a CGS
# # geneCGS = geneAll.copy()
# CM = mu.CMapMongo()
# CGSsigs = CM.find({'pert_type':'trt_sh.cgs','pert_iname':{'$in':list(geneAll)}},{'sig_id':True,'pert_iname':True},toDataFrame = True)

### single example
### use pert_explorer clustering tool
for cell in cellSet:
    out = wkdir + '/' + cell
    if not os.path.exists(out):
        os.mkdir(out)
    # examine drug rosiglitazone, on CPC006 plates, gold signatures
    # query = {'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True, 'pert_dose':{'$gt':3}}
    query = {'pert_type':'trt_sh.cgs','cell_id':cell,'pert_iname':{'$in':list(geneAll)}}

    self = PertExplorer(pert_id = None,
                        metric = 'wtcs',
                        query = query,
                        out = out)
    # cluster the data, store it on the PertExplorer instance
    self.hclust = HClust(self.score.copy(), self.annots.copy(), out = self.out,
                         symmfun = np.max, cutoff = 0.5)
    self.hclust.cluster()
    self.hclust.get_cluster_order()
    # self.hclust.draw_dendrogram(orientation = 'right')
    # self.hclust.save_dendrogram()
    clust_order = self.hclust.cluster_order.order().index
    # add cluster assignment to annotations
    self.annots['hclust_assignment'] = self.hclust.cluster_assignment
    # draw heatmap, ordered by the clustering; save
    #does showfig flag work?
    self.draw_heatmap(order = clust_order,
                      vmin = 0, vmax = 1,
                      cmap = cm.OrRd,
                      col_label = 'cell_id',
                      row_label = None,
                      annots_top = [('hclust_assignment', False)],
                      annots_right = ['pert_dose',
                                      ('distil_cc_q75', {'vmin' : 0, 'vmax' : 1}),
                                      ('distil_ss', {'vmin' : 0, 'vmax' : 10})],
                      show_title = False,
                      title = cell + ' Clustering for MTOR pathway',
                      showfig = False)
                      # annots_bottom = [('cell_lineage', True)],
    self.heatmap.save(format = 'png')
    # draw baseline data; save
    # self.showBaselines(order = self.annots.cell_id[clust_order].values,
    #                    clust_assignment = self.hclust.cluster_assignment.order().values)
    # self.saveBaselines(format = 'png')
    # write html file and annotations
    # self.write_html(extra_figs = [(self.hclust.dendro_fname, 'dendrogram')])
    self.write_annotations(encoding = 'latin-1')
    genesTested = set(self.annots['pert_iname'].values)
    genesNotTested = geneAll.difference(genesTested)
    print 'genes from this pathway not tested in ' + cell + ':'
    for x in genesNotTested:
        print x 
    skippedF = out + '/pathway_genes_skipped_' + cell +'.txt'
    gnt = pd.Series(list(genesNotTested))
    gnt.to_csv(skippedF, index=False, sep='\n')
    #### correlate random CGS pairs
    # retrieve random CGS sig ids from db
    CM = mu.CMapMongo()
    cellSigs = CM.find({'pert_type':'trt_sh.cgs','cell_id':cell},{'sig_id':True},toDataFrame = True)
    cellsigFrm = cellSigs['sig_id']
    n_cellSig = len(cellSigs)
    n_rands = 200
    # iRand = np.random.randint(0,n_cellSig,n_rands)
    iRand = np.random.choice(n_cellSig,n_rands,replace=False)
    high=None
    randSigsFrm = cellsigFrm.ix[iRand]
    randSigs = randSigsFrm.values
    queryer = Queryer()
    outDir = '/'.join([wkdir,'null_queries',cell])
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    queryer.set_params(out = outDir,
                        sig_id = randSigs.tolist(),
                        cid = randSigs.tolist(),
                        metric = self.metric,
                        row_space = 'lm')
    qr = queryer.run_query()
    randFrm = qr.score.loc[randSigs,randSigs]    
    # make array of random connections
    ilRand = np.triu_indices(len(randFrm),k=0)
    upRand = randFrm.values.copy()
    upRand[ilRand] = np.nan
    randFlat = upRand.flatten()
    uniqRand = upRand[~np.isnan(upRand)]
    # make array of observed connections
    il1 = np.triu_indices(len(self.score),k=0) #indices of upper triangle
    upScores = self.score.values.copy()
    upScores[il1] = np.nan 
    upFlat = upScores.flatten()
    uniqScores = upFlat[~np.isnan(upFlat)]
    ### make histogram
    plt.close()
    h1 = plt.hist(uniqRand,50,color='b',range=[-1,1],label=['random connections'],normed=True,alpha=0.5)
    h2 = plt.hist(uniqScores,50,color='r',range=[-1,1],label='within pathway connections',normed=True,alpha=0.5)
    plt.legend()
    plt.xlabel(self.metric)
    plt.ylabel('relative freq')
    plt.title(cell + ' distribution of pathway connections\n'+pathName)
    plt.savefig(os.path.join(out,cell + '_pathway_connect_dist.png'))
    plt.close()

#reload pert_explorer content
import cmap.analytics.pert_explorer as pertE
reload(pertE)
from cmap.analytics.pert_explorer import PertExplorer
