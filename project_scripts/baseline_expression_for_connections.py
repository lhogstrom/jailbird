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

### check the baseline expression of target genes in cell lines. Does this predict
# the connection of the drug?

### asumption: 
# if a cell line does not express the expected target
# then any expression consequences of the drug are considered
# off target


##########################################
### 1 check the ratio of is gold #########
### for cell lines expressing target #####
##########################################

#get instances to all drugs that have a target and are in summly space
qr = MC.sig_info.find({'pert_id':{'$in':po.cpTested}},
        {'sig_id':True,'pert_iname':True,'pert_id':True,'is_gold':True,'cell_id':True},
        toDataFrame=True)

goldDrugGrped = qr.groupby(['pert_id','is_gold'])
drugGrped = qr.groupby('pert_id')

#pick a gene find its expression
# set dummy frame to establish the index
goldFrame = pd.DataFrame(np.arange(6),columns=['dummy']) 
goldFrame = goldFrame.reindex([(False, False), (False, True), (True, False), (True, True), (u'not_measured', False), (u'not_measured', True)])
# goldFrame.index = [(False, False), (False, True), (True, False), (True, True), (u'not_measured', False), (u'not_measured', True)]
for gene in oneFrm.index:
    # gene='BCL2'
    drugs = summDict[gene]
    for drug in drugs: 
        # drug = drugs[0]
        exprDict = oneFrm.ix[gene].values[0]
        exprSer = pd.Series(oneFrm.ix[gene].values[0])
        exprSer.name = gene
        #make a frame for each drug
        iDrug = drugGrped.groups[drug]
        drugFrm = qr.ix[iDrug]
        drugFrm['baseline_expr'] = 'not_measured'
        #fill in expression status of target gene
        cellGrped = drugFrm.groupby('cell_id')
        for cell in cellGrped.groups.keys():
            icell = cellGrped.groups[cell]
            if exprDict.has_key(cell):
                if exprDict[cell]:
                    drugFrm.ix[icell,'baseline_expr'] = True
                else:
                    drugFrm.ix[icell,'baseline_expr'] = False
                # drugFrm.ix[icell,'baseline_expr'] = 'cheese'
        # baseLinGrped = drugFrm.groupby('baseline_expr')
        # count number of gold in expressing
        baseGoldGrped = drugFrm.groupby(['baseline_expr','is_gold'])
        sumSer = baseGoldGrped['is_gold'].size()
        sumSer.name = gene + ':' + drug
        sumFrm = pd.DataFrame(sumSer)
        goldFrame = pd.concat([goldFrame,sumFrm],axis=1)
        # if goldFrame.shape[1] < 2:
        #     goldFrame = goldFrame.reindex([(False, False), (False, True), (True, False), (True, True), (u'not_measured', False), (u'not_measured', True)])
goldFrame.__delitem__('dummy')
# convert nan to zeros
goldFrame = goldFrame.replace(np.nan,0)

# what is the percentage of is_gold in cell lines not expressing target gene
notExprNotGold = goldFrame.ix[0,:] #counts of not gold in cell lines where target is not expressed
notExprIsGold = goldFrame.ix[1,:] #counts of is_gold in cell lines where target is not expressed
notExprPercGold = notExprIsGold/ (notExprIsGold + notExprNotGold) * 100

# what is the percentage of is_gold in cell lines expressing target gene
ExprNotGold = goldFrame.ix[2,:] #counts of not gold in cell lines where target gene is expressed
ExprIsGold = goldFrame.ix[3,:] #counts of is_gold in cell lines where target is expressed
ExprPercGold = ExprIsGold/ (ExprIsGold + ExprNotGold) * 100

sortEPG = ExprPercGold.order(ascending=False)
sortNEPG = notExprPercGold.reindex(sortEPG.index)

plt.plot(sortEPG,'.r',label='target gene expressed')
plt.plot(sortNEPG,'.b',label='target gene not expressed')
plt.xlabel('drugs-target pairs')
plt.ylabel('percent is_gold signatures')
plt.legend(loc="upper right")
outF = '/xchip/cogs/hogstrom/analysis/scratch/baseline_of_drug_targets.png'
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()



#what if a gene has multiple targets? which one is most important?
