#! /usr/bin/env python
'''
examine nonsynonymous mutation list proveded by Aaron
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

import numpy as np
from cmap.analytics.cluster import HClust
import cmap
from os import path
from matplotlib import cm

wkdir = '/xchip/cogs/projects/oncoDome/mutation_clustering'

gl = ['APC', 
    'ATM', 
    'ATR', 
    'BRAF', 
    'CBL', 
    'CHEK2', 
    'EP300', 
    'FANCM', 
    'KRAS', 
    'NRAS', 
    'PIK3CA', 
    'SMARCA4', 
    'TP53', 
    'TP63']

#which of these drugs have genes targeting them?
#is clustering of KD itself 

CM = mu.CMapMongo()    
# CGSsigs = CM.find({'pert_type':'trt_sh.cgs','pert_iname':{'$in':list(geneAll)}},{'sig_id':True,'pert_id':True},toDataFrame = True)
CGSsigs = CM.find({'pert_type':'trt_sh.cgs','pert_iname':'KRAS'},{'sig_id':True,'pert_id':True},toDataFrame = True)

CM = mu.CMapMongo()
pert_id_dict = {}
for gene in gl:
    pID = CM.find({'pert_type':'trt_sh.cgs','pert_iname':gene},{'pert_id':True},limit=1)
    if pID:
        pert_id_dict[gene] = pID[0]

for gene in pert_id_dict:
    out = wkdir + '/' + gene
    if not os.path.exists(out):
        os.mkdir(out)
    query = {'pert_type':'trt_sh.cgs'} 
    self = PertExplorer(pert_id = pert_id_dict[gene],
                        metric = 'spearman',
                        query = query,
                        out = out)
    self.hclust = HClust(self.score.copy(), self.annots.copy(), out = self.out,
                         symmfun = np.max, cutoff = 0.5)
    self.hclust.cluster()
    self.hclust.get_cluster_order()
    self.hclust.draw_dendrogram(orientation = 'right', showfig=False)
    self.hclust.save_dendrogram()
    clust_order = self.hclust.cluster_order.order().index
    # add cluster assignment to annotations
    self.annots['hclust_assignment'] = self.hclust.cluster_assignment
    # draw heatmap, ordered by the clustering; save
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
                      showfig=False,
                      title = 'Clustering for the knockdown of ' + gene)
    self.heatmap.save(format = 'png')
    self.showBaselines(order = self.annots.cell_id[clust_order].values,
                       clust_assignment = self.hclust.cluster_assignment.order().values,
                       showfig=False)
    self.saveBaselines(format = 'png')
    # write html file and annotations
    self.write_html(extra_figs = [(self.hclust.dendro_fname, 'dendrogram')])
    self.write_annotations(encoding = 'latin-1')

# #initiate baseline object
# BaseLine.__init__(self, self.pert_id, self.out, 
#               cell_id = np.sort(self.annots.cell_id.unique()))




