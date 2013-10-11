'''
PertExplorer class to analyze avicin compounds 

July 2013
'''
import numpy as np
from cmap.analytics.pert_explorer import PertExplorer
from cmap.analytics.cluster import HClust
import cmap
import os
from os import path
from matplotlib import cm
import cmap.util.mongo_utils as mu
import subprocess
import cmap.tools.sig_dose_tool as sdt
import cmap.io.gct as gct

# get directory
wkdir = '/xchip/cogs/projects/avicins' 
out = wkdir + '/pert_explorer_avicins_Oct8'
if not os.path.exists(out):
    os.mkdir(out)

avicinsBrds = ['BRD-A15100685','BRD-A33746814','BRD-A69592287','BRD-A70150975'] #avicin-d, avicin-g, oxetane, hydroxyl,
pDescDict = {'BRD-A15100685':'avicin-d','BRD-A33746814':'avicin-g','BRD-A69592287':'oxetane','BRD-A70150975':'hydroxyl'}

for brd in avicinsBrds:
    # out = wkdir + '/pert_explorer_avicins/' + brd
    # if not os.path.exists(out):
    #     os.mkdir(out)
    # examine drug rosiglitazone, on CPC006 plates, gold signatures
    # query = {'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True} 
    query = {'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True, 'pert_dose':{'$gt':3}}
    self = PertExplorer(pert_id = brd,
                        metric = 'wtcs',
                        query = query,
                        out = out)
    # cluster the data, store it on the PertExplorer instance
    self.hclust = HClust(self.data.copy(), self.annots.copy(), out = self.out,
                         symmfun = np.max, cutoff = 0.5)
    self.hclust.cluster()
    self.hclust.get_cluster_order()
    self.hclust.draw_dendrogram(orientation = 'right')
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
                      showfig = False,
                      title = 'Example clustering for sig_id ' + pDescDict[brd])
    self.heatmap.save(format = 'png')
    # draw baseline data; save
    self.showBaselines(order = self.annots.cell_id[clust_order].values,
                       clust_assignment = self.hclust.cluster_assignment.order().values)
    self.saveBaselines(format = 'png')
    # write html file and annotations
    self.write_html(extra_figs = [(self.hclust.dendro_fname, 'dendrogram')])
    self.write_annotations(encoding = 'latin-1')

# make mongo query, write sig_ids to file
max_processes=7
processes = set()
for brd in avicinsBrds:
  CM = mu.CMapMongo()
  qRes = CM.find({'pert_id':{'$regex':brd},'pert_dose':{'$gt':3}},
              {'sig_id':True},
              toDataFrame=True)
  iname = pDescDict[brd]
  outF = out + '/' + iname + '.grp'
  qRes.to_csv(outF,index=False,header=False)
  #run sig_introspect
  cmd = ' '.join(['rum -q hour',
         '-d sulfur_io=100',
         '-o ' + out,
         '-x sig_introspect_tool ',
         '--sig_id ' + outF,
         '--metric wtcs',
         '--out ' + out])
  # os.system(cmd)
  # processes.add(subprocess.Popen(cmd,shell=True))
  # if len(processes) >= max_processes:
  #     os.wait()
  #     processes.difference_update(
  #         p for p in processes if p.poll() is not None)
  # make gct file for sig_dose
  afPath = cmap.score_path
  gt = gct.GCT()
  gt.read(src=afPath,cid=list(qRes['sig_id'].values),rid='lm_epsilon')
  outGctPref = out + '/' + iname + 'dose_data'
  gt.write(outGctPref)
  gt.matrix.shape
  outGct = outGctPref + '_n' + str(gt.matrix.shape[1]) + 'x' + '978.gctx' # what is the full file name
  ###sig dose tool
  cmdSDT = ' '.join(['python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/sig_dose_tool.py ' + outGct,
         '--out ' + out])
  processes.add(subprocess.Popen(cmdSDT,shell=True))
  if len(processes) >= max_processes:
      os.wait()
      processes.difference_update(
          p for p in processes if p.poll() is not None)

### find estrodiol instances
CM = mu.CMapMongo()
qRes = CM.find({'pert_iname':{'$regex':'stradiol'},'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True},
            {'sig_id':True,'pert_iname':True,'pert_id':True},
            toDataFrame=True)
brdSet = set(qRes['pert_id'].values)

query = {'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True, 'pert_dose':{'$gt':3}}
self = PertExplorer(pert_id = brd,
                    metric = 'wtcs',
                    query = query,
                    out = out)
