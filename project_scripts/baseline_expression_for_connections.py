'''
look up baseline expression of genes - does this help to 
predict connection in a given set of cell lines

November 2013
'''
import numpy as np
from cmap.analytics.pert_explorer import PertExplorer
import cmap.analytics.pcla as pcla
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
from scipy import stats
import matplotlib.pyplot as plt
import test_modules.load_TTD_drug_class as ldc

# get directory
# dir1 = '/xchip/cogs/projects/pharm_class' 
# wkdir = dir1 + '/estrodiol_analysis_Oct11'
# if not os.path.exists(wkdir):
#     os.mkdir(wkdir)

### load in data for individual groups
llo = ldc.label_loader()
pclDict = llo.load_drugbank_by_gene(group_by_action=False)

#get updated pcl dict (to summly space) using pcla module
### class analysis
wkdir = '/xchip/cogs/projects/pharm_class/TTd_Oct29'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
reload(pcla)
metric = 'wtcs'
po = pcla.PCLA(pclDict,    
                    metric,
                    wkdir,
                    pairwise_prefix='pairwise_matrices',
                    rankpt_metric='mean_rankpt_4',
                    sum_score_metric='sum_score_4',
                    row_space = 'lm',
                    cell_match_mode=True)
po.get_inames()
po.load_summly_mtrx()
summDict = po.pclResultDict
dbGenes = summDict.keys()
dbGenes.remove('-666')

#look up baseline expression of these genes
MC = mu.MongoContainer()
ci = MC.gene_info.find({'pr_gene_symbol':{'$in':dbGenes}},
            {'is_expressed':True,'pr_gene_symbol':True},
            toDataFrame=True)
gsGrped = ci.groupby('pr_gene_symbol')
oneFrm = gsGrped.first() #remove dup data if a gene has more than one affy probe

# define core cell line
cellList = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
#find which genes are expressed in about 
#half of the cell line
halfExpressingFrm = pd.DataFrame() 
for gene in oneFrm.index:
    exprSer = pd.Series(oneFrm.ix[gene].values[0])
    coreSer = exprSer.ix[cellList]
    coreSer.name = gene
    if (coreSer.sum() >= 3) and (coreSer.sum() <= 7):
        coreFrm = pd.DataFrame(coreSer).T
        targetDrugs = summDict[gene]
        coreFrm['drugs_targeting'] = ':'.join(targetDrugs)
        halfExpressingFrm = pd.concat([halfExpressingFrm,coreFrm],axis=0)
outF = '/xchip/cogs/hogstrom/analysis/scratch/baseline_of_drug_targets_in_core_cell_lines.txt'
halfExpressingFrm.to_csv(outF,headers=True,sep='\t')