'''
Find data related to doxyrubicin

5/15/2014
'''
import pandas as pd
import os
import glob
import shutil
import matplotlib.pyplot as plt
import numpy as np

import cmap.analytics.NMF_benchmarks as nmfb
import cmap.util.mongo_utils as mu
import cmap.io.gmt as gmt
import cmap.io.gct as gct

wkdir = '/xchip/cogs/hogstrom/bathe/gordonov/doxorubicin'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

### get doxyrubicin signatures
mc = mu.MongoContainer()
## KD
doxFrm = mc.sig_info.find({'pert_iname':'doxorubicin','pert_type':'trt_cp'},
            {'dn100_full':False,'up100_full':False,'dn100_bing':False,'up100_bing':False,'dn50_lm':False,'up50_lm':False,'pert_idose':False,'pert_dose_unit':False},toDataFrame=True)
## OE
# doxFrm = mc.sig_info.find({'pert_iname':'doxorubicin','pert_type':'trt_cp'},
#             {'dn100_full':False,'up100_full':False,'dn100_bing':False,'up100_bing':False,'dn50_lm':False,'up50_lm':False,'pert_idose':False,'pert_dose_unit':False},toDataFrame=True)
cellGrped = doxFrm.groupby('cell_id')

# dose list for each cell line
cellGrped.apply(lambda x: set(x.pert_dose))

######################
### run query tool ###
######################

sigs = doxFrm.sig_id
sig_file = os.path.join(wkdir,'doxorubicin_sig_ids.grp')
sigs.to_csv(sig_file,index=None)
metric = 'wtcs'
queryDir = os.path.join(wkdir,'sig_query')
if not os.path.exists(queryDir):
    os.mkdir(queryDir)
cmd = ' '.join(['rum -q local -f sig_query_tool',
         '--sig_id ' + sig_file,
         '--metric ' + metric,
         '--column_space full',
         '--out ' + queryDir,
         '--mkdir false',             
         '--save_tail false'])
# os.system(cmd)

######################
### load query results ###
######################

query_res = os.path.join(queryDir,'result_WTCS.LM.COMBINED_n150x476251.gctx')
gt = gct.GCT(query_res)
gt.read()
qRes = gt.frame

resSer = qRes['CPC004_A375_6H:BRD-A52530684-001-01-1:10']
resSer.sort()

######################
### plot dox signature counts ###
######################

### dox signature counts
cellGrped = doxFrm.groupby('cell_id')
counts = cellGrped.cell_id.count()
# counts.plot(kind='bar',x='cell_line', y='number_of_signatures')
counts.plot(kind='bar',title='doxorubicin signature counts')
outF = os.path.join(wkdir,'dox_signature_counts.png')
plt.savefig(outF)
plt.close()

######################
### TWIST OE sigs? ###
######################

emtTFs = ['TWIST1', 'SNAI1', 'SNAI2']
tfFrm = mc.sig_info.find({'pert_iname':{'$in':emtTFs},'pert_type':'trt_oe'},
            {'dn100_full':False,'up100_full':False,'dn100_bing':False,'up100_bing':False,'dn50_lm':False,'up50_lm':False,'pert_idose':False,'pert_dose_unit':False},toDataFrame=True)
tfFrm.ix[:,['pert_iname','cell_id','distil_cc_q75','pert_itime']]


