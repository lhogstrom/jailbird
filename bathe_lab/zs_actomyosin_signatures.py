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

wkdir = '/xchip/cogs/hogstrom/bathe/gordonov'
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

# calculate pairwise correlations
corrMtrx = np.corrcoef(ds,rowvar=0)
corrFrm = pd.DataFrame(corrMtrx, index=ds.columns,columns=ds.columns)

graphDir = os.path.join(wkdir,gene_list_name + '_OE_replicates')
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
pairwiseDict = {}
cell_grped = cgsFrm.groupby(['pert_iname','cell_id'])
cell_lines = set(cgsFrm.cell_id)
for cell in cell_lines:
    for grp in cell_grped:
        if grp[0][1] == cell:
            frm = grp[1]
            frm.shape
            id_list = [item for sublist in frm.distil_id for item in sublist]
            mtch_corr = corrFrm.reindex(index=id_list,columns=id_list)
            il = np.triu_indices(len(mtch_corr),k=0)
            upMeans = mtch_corr.values.copy()
            upMeans[il] = np.nan
            randFlat = upMeans.flatten()
            uniqRand = upMeans[~np.isnan(upMeans)]
            pairwiseDict[grp[0][0]] = uniqRand
    # sort by median
    inSer = pd.Series(pairwiseDict)
    inMedian = inSer.apply(np.median)
    inMedian.sort()
    inMedian = inMedian[~np.isnan(inMedian)]
    inSorted = inSer[inMedian.index]
    ### make boxplot
    fig = plt.figure(figsize=(8, 14), dpi=50)
    plt.boxplot(inSorted,vert=0)
    plt.xlim((-1,1))
    tickList = [x for x in inSorted.index]
    plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
    plt.tick_params(labelsize=8)
    plt.xlabel('Pearson corr',fontweight='bold')
    plt.title(cell + ' - pairwise hairpin correlations',fontweight='bold')
    outF = os.path.join(graphDir,cell+'_pairwise_comparison_boxplot_cytoskeleton_OE.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close() 

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
    cellFrm = ds_lung.ix[:,sigs.values]
    nGt = gct.GCT()
    nGt.build_from_DataFrame(cellFrm)
    outF = cellDir + '/' + cell + '_TA_JUN10_' + processesed_type
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



# create ranked lists (connections)
# rnkpt, mutual information



