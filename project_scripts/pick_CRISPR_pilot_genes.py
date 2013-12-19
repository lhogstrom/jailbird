'''
decide on subset of genes to pilot with CRISPRS (~360)
'''
import cmap.analytics.pcla as pcla
import numpy as np
import os
import cmap.io.gmt as gmt
import cmap
import pandas as pd
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu

# Overall criteria for selection - Include the following (overlap is ok):
# 1) 75-100 LM
# 2) 50-100 genes targeted by small molecule
# 3) 25-50 genes with 6+ shRNAs
# 4) 50 not expressed in some LINCS core line cells
# 5) 50 with really good shRNAs - good KD and good hairpin consensus
# 6) 50 with bad shRNAs -no KD and bad hairpin consensus

wkdir = '/xchip/cogs/projects/CRISPR/pilot'
###########################
### target KD efficiency #
###########################

#get differential expression of the shRNAs
#get all CGS that target a LM
mc = mu.MongoContainer()
cgsInfo = mc.sig_info.find({'target_is_lm':True,'pert_type':'trt_sh.cgs'},
            {},toDataFrame=True)
geneGrped = cgsInfo.groupby('pert_iname')
medianTargetExpr = geneGrped['target_zs'].median()
medianTargetExpr.sort()

plt.hist(medianTargetExpr)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('target_zs (median across all cell lines)',fontweight='bold')
plt.title('LM genes targetd by shRNA - differential expression of target')
outF = os.path.join(wkdir, 'target_expression_LM.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

KDwell = medianTargetExpr[medianTargetExpr < -5]
KDbad = medianTargetExpr[medianTargetExpr > 0]

###########################
### shRNA counts ##########
###########################

cgsFrm = mc.sig_info.find({'pert_type':'trt_sh.cgs'},
            {'dn100_full':False,'up100_full':False,'dn100_bing':False,'up100_bing':False,'dn50_lm':False,'up50_lm':False},toDataFrame=True)
geneGrped= cgsFrm.groupby('pert_iname')
medNsample = geneGrped['distil_nsample'].median()
medNsample.sort(ascending=False)
largeNsample = medNsample[medNsample>10]

###########################
### baseline expr. ########
###########################

#baseline expression in the core cell lines
#dynamic expression range

ci = mc.gene_info.find({'is_lm':True,'pr_pool_id':'epsilon'},
            {'is_expressed':True,'pr_gene_symbol':True},
            toDataFrame=True)
gsGrped = ci.groupby('pr_gene_symbol')
oneFrm = gsGrped.first()

###########################
### targets with Drugs ####
###########################

#drug targets - important ones, well connecting ones --> plus all LM
filePCLgrps = '/xchip/cogs/projects/pharm_class/pcl_shared_target.txt'
pclFrm = pd.io.parsers.read_csv(filePCLgrps,sep='\t')
drugFrm = pclFrm[pclFrm['src'] == 'DRUG_BANK']
targetList = drugFrm['class'].values
targetSet = set(targetList)
#how many of these have a cgs in CMAP
targetFrm = mc.sig_info.find({'pert_type':'trt_sh.cgs','pert_iname':{'$in':list(targetSet)}},
            {'sig_id':True,'pert_iname':True},toDataFrame=True)
cgsTargetdSet = set(targetFrm['pert_iname'])


#crisprFrm 

# include these guys:
largeNsample.index #could add more of these
cgsTargetdSet
KDwell
KDbad

#get all lm genes 
# mc = mu.MongoContainer()
# geneInfo = mc.gene_info.find({'is_lm':True,'pr_pool_id':'epsilon'},
#             {'pr_gene_symbol':True,'pr_id':True,'is_l1000':True},toDataFrame=True)
# #check that it doesn't have a known pert_iname
# symbolSer = pd.Series(geneInfo['pr_gene_symbol'],index=geneInfo['pr_id'])
# symbolSer = pd.Series(geneInfo['pr_gene_symbol'])
# symbolSer.index=geneInfo['pr_id']
