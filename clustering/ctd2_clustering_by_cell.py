#! /usr/bin/env python
'''
take every ctd2 signature from one cell line
put them in a big matrix and cluster them
repeat for each cell line

which compounds stay together most closely and consistantly? - Average distance of drug-drug pairs
does this reflect class labels?

histogram:
cluster DMSOs Is the average pairwise correlation higher or lower?
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

wkdir = '/xchip/cogs/projects/target_id/ctd2_cp_clustering'
if not os.path.exists(wkdir):
    os.makedirs(wkdir)

#what are the cps that were tested in all 50 cell lines?
CM = mu.CMapMongo()    
qr = CM.find({'pert_type':'trt_cp','sig_id' : {'$regex' : 'CPC006'}},{'cell_id':True,'pert_id':True},toDataFrame = True)
ctdCells = set(qr['cell_id'])
#cps from non-core cell line - LOVO
cpQ = CM.find({'pert_type':'trt_cp','cell_id':'LOVO','sig_id' : {'$regex' : 'CPC006'}},{'cell_id':True,'pert_id':True},toDataFrame = True)
ctdCps = list(cpQ['pert_id'])

metric = 'spearman'
cell = 'LOVO'

out = wkdir + '/' + cell
if not os.path.exists(out):
    os.mkdir(out)
query = {'pert_type':'trt_cp','cell_id':cell,
        'pert_id':{'$in':ctdCps}}
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
PE.hclust.draw_dendrogram(orientation = 'right')
PE.hclust.save_dendrogram()
### add cluster assignment to annotations
PE.annots['hclust_assignment'] = PE.hclust.cluster_assignment
### draw heatmap, ordered by the clustering; save
PE.draw_heatmap(order = clust_order,
                  vmin = 0, vmax = 1,
                  cmap = cm.OrRd,
                  col_label = 'pert_iname',
                  row_label = 'pert_iname',
                  annots_top = [('hclust_assignment', False)],
                  annots_right = ['pert_dose',
                                  ('distil_cc_q75', {'vmin' : 0, 'vmax' : 1}),
                                  ('distil_ss', {'vmin' : 0, 'vmax' : 10})],
                  show_title = False,
                  title = cell + ' Clustering ctd2',
                  ticklabel_fontsize = 4,
                  showfig = False)
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



# genesTested = set(PE.annots['pert_iname'].values)
# geneAll = set(self.inputGeneSet)
# genesNotTested = geneAll.difference(genesTested)
# print 'genes from this pathway not tested in ' + cell + ':'
# for x in genesNotTested:
#     print x 
# skippedF = out + '/pathway_genes_skipped_' + cell +'.txt'
# gnt = pd.Series(list(genesNotTested))
# gnt.to_csv(skippedF, index=False, sep='\n')

PE.heatmap.row_labels






