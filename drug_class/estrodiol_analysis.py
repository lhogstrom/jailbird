'''
analyze estrodiol compounds 

October 2013
'''
import numpy as np
from cmap.analytics.pert_explorer import PertExplorer
from cmap.analytics.cluster import HClust
import cmap.analytics.sc as sc
import cmap
import os
from os import path
from matplotlib import cm
import cmap.util.mongo_utils as mu
import subprocess
import cmap.tools.sig_dose_tool as sdt
import cmap.io.gct as gct
import pandas as pd
import cmap.io.gmt as gmt

# get directory
dir1 = '/xchip/cogs/projects/pharm_class' 
wkdir = dir1 + '/estrodiol_analysis_Oct11'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

pDescDict = {'BRD-A60070924':'17-alpha-estradiol',
    'BRD-K04218075':'clomiphene',
    'BRD-K95309561':'dienestrol',
    'BRD-K45330754':'diethylstilbestrol',
    'BRD-K04046242':'equilin',
    'BRD-A74907996':'equol',
    'BRD-K18910433':'estradiol-17-beta',
    'BRD-K17016787':'estriol',
    'BRD-K81839095':'estrone',
    'BRD-A83237092':'fulvestrant',
    'BRD-K43797669':'genistein',
    'BRD-K63828191':'raloxifene',
    'BRD-K93754473':'tamoxifen',
    'BRD-K67174588':'toremifene'}
estrodiolBrds = pDescDict.keys()

#
for brd in estrodiolBrds:
    out = wkdir + '/' pDescDict[brd]
    if not os.path.exists(out):
        os.mkdir(out)
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

### run sig_introspect on each compound seperatly
max_processes=7
processes = set()
sigIdDict = {}
# make S-C plots, color by cell line
for brd in estrodiolBrds:
    # out = wkdir + '/' + pDescDict[brd]
    # if not os.path.exists(out):
    #     os.mkdir(out)
    CM = mu.CMapMongo()
    #all cell lines
    # qRes = CM.find({'pert_id':brd,'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True},
    #         {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
    #         toDataFrame=True)
    #MCF7 only query
    qRes = CM.find({'pert_id':brd,'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True,'cell_id':'MCF7'},
            {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
            toDataFrame=True)
    if type(qRes) == int:
        continue
    sigIDs = list(qRes['sig_id'].values)
    sigIdDict[brd] = sigIDs
    # make SC plots
    out = wkdir + '/SC_plots_MCF7'
    if not os.path.exists(out):
        os.mkdir(out)
    sco = sc.SC()
    sco.add_sc_from_mongo_sig_ids(sigIDs)
    scOut = out + '/' + pDescDict[brd] + '_sc.png'
    sco.plot(out=scOut)
    ### run sig_intorspect
    outIntrospect = out + '/brew_introspect_results_MCF7'
    # outIntrospect = out + '/brew_introspect_results_all_cell_lines'
    if not os.path.exists(outIntrospect):
        os.mkdir(outIntrospect)
    qSer = qRes['sig_id']
    iname = pDescDict[brd]
    outF = outIntrospect + '/' + iname + '_sig_ids.grp'
    qSer.to_csv(outF,index=False,header=False)
    #run sig_introspect
    cmd = ' '.join(['rum -q hour',
         '-d sulfur_io=100',
         '-o ' + outIntrospect,
         '-x sig_introspect_tool ',
         '--sig_id ' + outF,
         '--metric wtcs',
         '--out ' + outIntrospect])
    os.system(cmd)
    # processes.add(subprocess.Popen(cmd,shell=True))
    # if len(processes) >= max_processes:
    #     os.wait()
    #     processes.difference_update(
    #         p for p in processes if p.poll() is not None)

# sim_set - pert_mfc
cellList = ['HEPG2', 'A375', 'PC3', 'MCF7', 'HA1E', 'HT29', 'HCC515', 'A549']
for cell in cellList:
    CM = mu.CMapMongo()
    qRes = CM.find({'pert_id':{'$in': estrodiolBrds},'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True,'cell_id':cell},
            {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
            toDataFrame=True)
    # look at pert_id-pert_mfc_id pairs
    pairs = qRes['pert_id'] + ':' + qRes['pert_mfc_id']
    pairSet = set(pairs.values)
    ### run sig_introspect on all of the sig_ids at once
    outIntrospect = wkdir + '/brew_introspect_all_cps_' + cell
    # outIntrospect = out + '/brew_introspect_results_all_cell_lines'
    if not os.path.exists(outIntrospect):
        os.mkdir(outIntrospect)
    qSer = qRes['sig_id']
    outF = outIntrospect + '/all_cps_' + cell + '_sig_ids.grp'
    qSer.to_csv(outF,index=False,header=False)
    #run sig_introspect
    cmd = ' '.join(['rum -q hour',
         '-d sulfur_io=100',
         '-o ' + outIntrospect,
         '-x sig_introspect_tool ',
         '--sig_id ' + outF,
         '--metric wtcs',
         '--out ' + outIntrospect])
    os.system(cmd)


# get instances from  sig_info - distil_id
CM = mu.CMapMongo()
sigList = qRes['sig_id']
sigInfo = mu.sig_info(sigList)
distilDict = {}
for sig in sigInfo:
    distilDict[sig['sig_id']] = sig['distil_id']
# # loop through each compound and write grp of roast ids
# for brd in estrodiolBrds:
#     iname = pDescDict[brd]
#     # cp dir
#     out = wkdir + '/' + iname
#     if not os.path.exists(out):
#         os.mkdir(out)
#     # introspect dir
#     outIntrospect = out + '/roast_introspect_MCF7'
#     if not os.path.exists(outIntrospect):
#         os.mkdir(outIntrospect)        
#     if brd in sigIdDict:
#         sigs = sigIdDict[brd]
#         distilIDList = []
#         for sig in sigs:
#             distilIDList.extend(distilDict[sig])
#         dIDSer = pd.Series(distilIDList)
#         outF = outIntrospect + '/' + iname + '_distil_ids.grp'
#         dIDSer.to_csv(outF,index=False,header=False)
#         #run sig_introspect - with roast data
#         cmd = ' '.join(['rum -q hour',
#              '-d sulfur_io=100',
#              '-o ' + outIntrospect,
#              '-x sig_introspect_tool',
#              '--build_id custom',
#              '--score ' + cmap.roast_path,
#              '--cid ' + outF,
#              '--metric wtcs',
#              '--out ' + outIntrospect])
#         os.system(cmd)

### put distil_ids into gmt form --> each row is a different (sig_id? or pert_csf_id?)

### make list of dictionaries for making gmt file
CM = mu.CMapMongo()
sigIdDict = {}
for brd in estrodiolBrds:
    #all cell lines
    qRes = CM.find({'pert_id':brd,'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True},
            {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
            toDataFrame=True)
    #MCF7 only query
    # qRes = CM.find({'pert_id':brd,'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True,'cell_id':'MCF7'},
    #         {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
    #         toDataFrame=True)
    if type(qRes) == int:
        continue
    sigIDs = list(qRes['sig_id'].values)
    sigIdDict[brd] = sigIDs
gmtList = []
for brd in sigIdDict:
    print brd
    sigs = sigIdDict[brd]
    distilIDList = []
    for sig in sigs:
        distilIDList.extend(distilDict[sig])
    # dIDSer = pd.Series(distilIDList)
    gmtDict = {}
    gmtDict['id'] = brd
    gmtDict['desc'] = pDescDict[brd]
    gmtDict['sig'] = distilIDList
    gmtList.append(gmtDict)
gmtOut = wkdir + '/roast_ids_all_cell_lines.gmt'
gmt.write(gmtList,gmtOut)



#cell_info
# look up ESR1 and ESR2 baseline expression in each cell line
