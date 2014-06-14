'''
Find data related to Simon's thesis project

5/15/2014
'''
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.io.gct as gct
import glob
import shutil
import cmap.util.mongo_utils as mu
import matplotlib.pyplot as plt
import numpy as np
import subprocess

wkdir = '/xchip/cogs/hogstrom/bathe/gordonov/KD_NMF_Jun10'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
gene_list_name = 'kegg_actin_cytoskeleton'
# gene_list_name = 'apriori_actomyosin'

# load kegg pathways
file_kegg = '/xchip/cogs/hogstrom/bathe/gordonov/c2.cp.kegg.v4.0.symbols.gmt'
gt = gmt.read(file_kegg)
keggFrm = pd.DataFrame(gt)
GeneList = keggFrm[keggFrm.id == 'KEGG_REGULATION_OF_ACTIN_CYTOSKELETON'].sig.values
GeneList = list(GeneList[0])

###
aprioriList = ['RAC1',
'CDC42', 
'RHOA',
'ROCK1',
'RICS',
'RHOA',
'PRKCA',
'PIK3CA',
'ARPC1A',
'MAPK',
'ERK',
'MAPK14',
'CAPN4',
'CAPN1',
'CAPN2',
'PTK2',
'SRC',
'NgR1',
'LINGO1',
'p75',
'TROY',
'MYH3',
'MYH6',
'MYH7',
'MYH9',
'MYH11',
'MYO1A',
'MYO5A',
'MYO6 ',
'MYO7A',
'MYO15A']
GeneList.extend(aprioriList)

### genomic perturbation
# shRNA
mc = mu.MongoContainer()
## KD
cgsFrm = mc.sig_info.find({'pert_type':'trt_sh','pert_iname':{'$in':list(GeneList)}},
            {'dn100_full':False,'up100_full':False,'dn100_bing':False,'up100_bing':False,'dn50_lm':False,'up50_lm':False,'pert_idose':False,'pert_dose_unit':False},toDataFrame=True)
## OE
# cgsFrm = mc.sig_info.find({'pert_type':'trt_oe','pert_iname':{'$in':list(GeneList)}},
#             {'dn100_full':False,'up100_full':False,'dn100_bing':False,'up100_bing':False,'dn50_lm':False,'up50_lm':False,'pert_idose':False,'pert_dose_unit':False},toDataFrame=True)
cgsFrm.index = cgsFrm.sig_id
iname_grp = cgsFrm.groupby('pert_iname')
idLst = [item for sublist in cgsFrm.distil_id for item in sublist]

# load raw replicate signatures
# sig1 = 'KDC006_A549_96H_X1_B6_DUO52HI53LO:O04'
file_zs = '/xchip/cogs/data/build/a2y13q1/zspc_n1328098x22268.gctx'
gt = gct.GCT()
gt.read(src=file_zs,cid=idLst,rid='lm_epsilon')
ds = gt.frame
processesed_type = 'ZSPC_LM'

## expand nested distil_ids
rowDict = {}
for x in cgsFrm.iterrows():
    rowSer = x[1]
    d_ids = rowSer.distil_id
    for d_id in d_ids:
        rowDict[d_id] = rowSer
distilAnnt = pd.DataFrame(rowDict)
distilAnnt = distilAnnt.T

# calculate pairwise correlations
# corrMtrx = np.corrcoef(ds,rowvar=0)
# corrFrm = pd.DataFrame(corrMtrx, index=ds.columns,columns=ds.columns)
# graphDir = os.path.join(wkdir,gene_list_name + '_OE_replicates')
# if not os.path.exists(graphDir):
#     os.mkdir(graphDir)

## split up columns for counting
colSer = pd.Series(ds.columns)
colSer.name = 'sig_id'
colSplit = colSer.str.split('_')
# annotate acording to sig_id fields
colFrame = pd.DataFrame(colSer)
colFrame['plate'] = colSplit.apply(lambda x: x[0])
colFrame['cell_line'] = colSplit.apply(lambda x: x[1])
colFrame['tp'] = colSplit.apply(lambda x: x[2])
colFrame['rep'] = colSplit.apply(lambda x: x[3])
colFrame['well'] = colSplit.apply(lambda x: x[4])

################################
### make cell line gct files ###
################################

# save matrix for each cell line in OE experiments
cell_grped = colFrame.groupby('cell_line')
for grpT in cell_grped:
    cell = grpT[0]
    cellDir = wkdir + '/' + cell
    if not os.path.exists(cellDir):
        os.mkdir(cellDir)
    grp = grpT[1]
    sigs = grp.sig_id
    cellFrm = ds.ix[:,sigs.values]
    nGt = gct.GCT()
    nGt.build_from_DataFrame(cellFrm)
    outF = cellDir + '/' + cell + '_actomyosin_' + processesed_type
    nGt.write(outF,mode='gctx')

# convert gctx to gct
### run 'use Java-1.7' before running script ###
for cell in cell_grped.groups.keys():
    print(cell)
    cellDir = wkdir + '/' + cell
    outGCT = cellDir + '/' + cell
    globRes = glob.glob(outGCT+'*.gctx')
    print(globRes[0])
    cmd2 = 'convert-dataset -i ' + globRes[0]
    os.system(cmd2)

################################
### make gene signature gmt ###
################################

cellGrped = distilAnnt.groupby('cell_id')
for cellTup in cellGrped:
    cell = cellTup[0]
    cellFrm = cellTup[1]
    cellDir = wkdir + '/' + cell
    outF = cellDir + '/KD_annotations.txt'
    # reformat sig_id
    cellFrm['distil_id_original'] = cellFrm.index
    cellFrm['mod_sig_id'] = cellFrm.distil_id_original.str.replace(':','.')
    cellFrm.index = cellFrm.mod_sig_id
    cellFrm.to_csv(outF,sep='\t')
    ### make gene signature groups - gmt file
    mtch_field = 'pert_iname'
    geneGrped = cellFrm.groupby(mtch_field)
    gmtList = []
    for grp in geneGrped:
        gmtDictUp = {}
        gmtDictUp['id'] = grp[0]
        gmtDictUp['desc'] = grp[0]
        # gmtDictUp['desc'] = str(list(set(grp[1][mtch_field])))
        gmtDictUp['sig'] = list(grp[1].index.values)
        gmtList.append(gmtDictUp)
    gmtOut = cellDir + '/actomyosin_kd_distil_id.gmt'
    gmt.write(gmtList,gmtOut)

#########################
### Run NMF projection ##
#########################

# COMPZ.MODZ_SCORE
nComponents = 20
# dimDict = {}
# for grp in cell_grped:
#     dimDict[grp[0]] = 'n'+str(grp[1].shape[0])+'x978'
dimDict = {'A375': 'n1684x978',
 'A549': 'n1410x978',
 'ASC': 'n260x978',
 'HA1E': 'n1445x978',
 'HCC515': 'n1163x978',
 'HEK293T': 'n39x978',
 'HEKTE': 'n222x978',
 'HEPG2': 'n1287x978',
 'HT29': 'n1679x978',
 'JURKAT': 'n15x978',
 'MCF7': 'n2042x978',
 'NPC': 'n312x978',
 'PC3': 'n2292x978',
 'SHSY5Y': 'n30x978',
 'SKL': 'n48x978',
 'SW480': 'n168x978',
 'U2OS': 'n112x978',
 'VCAP': 'n2883x978'}

#specifications for subprocess
processes = set()
max_processes = 9 
### run jobs
for cell in cell_grped.groups.keys():
    print cell
    dim = dimDict[cell]
    arg1 = wkdir + '/' + cell # working directory
    arg2 = cell + '_actomyosin_' + processesed_type + '_' + dim
    cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/NMF_code_no_viz.v2.R', # 
         arg1,
         arg2])
    # os.system(cmd)
    processes.add(subprocess.Popen(cmd,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes if p.poll() is not None)







# create ranked lists (connections)
# rnkpt, mutual information



