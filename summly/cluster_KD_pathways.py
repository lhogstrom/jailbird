#! /usr/bin/env python
'make a roc curve to look at drug-gene connnections from summly results'

import numpy as np
import os
import cmap.util.mongo_utils as mu
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
from cmap.analytics.pert_explorer import PertExplorer
from cmap.analytics.cluster import HClust
import cmap.io.gmt as gmt

wkdir = '/xchip/cogs/projects/target_id/KD_pathway_clustering'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

#pathway annotations from reactome 
pathGMT = '/xchip/cogs/projects/target_id/KD_pathway_clustering/ReactomePathways.gmt'
pathwayDict = gmt.read(pathGMT)
mtorPath = pathwayDict[1398]
mtorGenes = mtorPath['sig']


#load in onco gene list
inFile = '/xchip/cogs/projects/oncoDome/OncoDome_genes.txt'
cgsList = []
with open(inFile,'rt') as f:
    for string in f:
        splt = string[:-1]
        cgsList.append(splt)
geneAll = set(cgsList)
# check to see if gene has a CGS
# geneCGS = geneAll.copy()
CM = mu.CMapMongo()
CGSsigs = CM.find({'pert_type':'trt_sh.cgs','pert_iname':{'$in':list(geneAll)}},{'sig_id':True,'pert_iname':True},toDataFrame = True)

### single example
### use pert_explorer clustering tool
out = wkdir + '/pert_explorer_onco/'
if not os.path.exists(out):
    os.mkdir(out)
# examine drug rosiglitazone, on CPC006 plates, gold signatures
# query = {'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True, 'pert_dose':{'$gt':3}}
# query = {'pert_type':'trt_sh.cgs','pert_iname':{'$in':list(geneAll)}}
self = PertExplorer(pert_id = 'CGS001-238',
                    metric = 'wtcs',
                    query = query,
                    out = out)
# self = PertExplorer(metric = 'wtcs',
#                     query = query,
#                     out = out)
# cluster the data, store it on the PertExplorer instance
self.hclust = HClust(self.score.copy(), self.annots.copy(), out = self.out,
                     symmfun = np.max, cutoff = 0.5)
self.hclust.cluster()
self.hclust.get_cluster_order()
self.hclust.draw_dendrogram(orientation = 'right')
self.hclust.save_dendrogram()
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
                  annots_bottom = [('cell_lineage', True)],
                  title = 'Example clustering for sig_id BRD-K06854232')
self.heatmap.save(format = 'png')
# draw baseline data; save
self.showBaselines(order = self.annots.cell_id[clust_order].values,
                   clust_assignment = self.hclust.cluster_assignment.order().values)
self.saveBaselines(format = 'png')
# write html file and annotations
self.write_html(extra_figs = [(self.hclust.dendro_fname, 'dendrogram')])
self.write_annotations(encoding = 'latin-1')