import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd

work_dir = '/xchip/cogs/hogstrom/analysis/summly/cp_class'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

drugFile = '/xchip/cogs/projects/target_id/ctd2_annots/ctd2_merged_mapped_genes.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
targetSet = set(drugLabels['gene_dw_annot'])
bigGroups = ['AKT','BRAF','EGFR','HDAC','MAPK','MTOR','HSP90','PI3K','PPARG','TP53']

brdGrpList = []
for grp in bigGroups:
    grpPerts = drugLabels['pert_id'][drugLabels['gene_dw_annot'] == grp]
    brdGrpList.extend(grpPerts.values)

allCtd2 = drugLabels['pert_id']

### run queries for cps
cgsCells = ['A375', 'A549', 'ASC', 'HA1E', 'HEPG2', 'HCC515', 'HT29', 'MCF7', 'NPC', 'PC3', 'VCAP']
processes = set()
file1 = work_dir + '/' + 'ctd2_sigs.grp'
for pert in allCtd2:
    sigs = []
    for cell in cgsCells:
        CM = mu.CMapMongo()
        #dose min
        # pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell,'pert_dose':{'$gt':5}},{'sig_id':True},limit=1)
        #no dose min
        pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell},{'sig_id':True},limit=1)
        if pert_query:
            sigs.append(pert_query[0])
        # else:
        #     print 'no sig for ' + pert + ' in ' + cell
    with open(file1, 'a') as f:
        [f.write(x + '\n') for x in sigs]
metric = 'wtcs'
queryDir = work_dir + '/sig_query_no_7371'
if not os.path.exists(queryDir):
    os.mkdir(queryDir)
cmd = ' '.join(['rum -q local -f sig_query_tool',
         '--sig_id ' + file1,
         '--metric ' + metric,
         '--column_space full',
         '--out ' + queryDir,
         '--mkdir false',
         '--save_tail false'])
         # '--row_space bing', 
os.system(cmd)


metric = 'wtcs'
queryDir = work_dir + '/test_sig_query'
if not os.path.exists(queryDir):
    os.mkdir(queryDir)
file1 = '/xchip/cogs/hogstrom/analysis/summly/cp_class/sigs_combo.grp'
cmd = ' '.join(['rum -q local -f sig_query_tool',
         '--sig_id ' + file1,
         '--metric ' + metric,
         '--column_space full',
         '--out ' + queryDir,
         '--mkdir false',
         '--save_tail false'])

