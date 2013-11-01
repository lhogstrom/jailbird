import cmap.analytics.gene_mod as gm 
import cmap
import test_modules.load_TTD_drug_class as ldc
import cmap.util.mongo_utils as mu
import pandas as pd
import os

# get directory
dir1 = '/xchip/cogs/projects/TRIB1' 
wkdir = dir1 + '/TRIB1_analysis_Oct21'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

## 
# gmo = gm.GeneMod()
# gmo.load_data_from_gctx(src=cmap.score_path,symbols='TRIB1')
# gmo.z_score_filter(4)
# gmo.signature_strength_filter(8)
# gmo.cell_type_filter('HEPG2')
# # gmo.sc_plot(out=wkdir+'sc_hits.png')
# gmo.scatter()
# gmo.expression_histogram()

# gmo.ssr_plot(out=wkdir+'ssr_hits.png')
# cid = [gmo.cid[x] for x in gmo.reg_ind]

#load in drugbank annotations
reload(ldc)
llo = ldc.label_loader()
pclDict = llo.load_TTD()
dbDict = llo.load_drugbank_by_gene(group_by_action=False)
geneTargets = dbDict.keys()
#put drug-gene relationships in a file
tupList = []
for gene in dbDict:
    # make tuple
    cps = dbDict[gene]
    for cp in cps:
        tup = (cp, gene)
        tupList.append(tup)
dbFrm = pd.DataFrame(tupList,columns=['brd','target_gene'])
dbSer = pd.Series(dbFrm['target_gene'])
dbSer.index = dbFrm['brd']
cpToGenDict = dbSer.to_dict()


#which drugs do we have that target lm or bing genes
test_space = 'is_lm'
test_space = 'is_bing'
cpSet = set(dbFrm['brd'])
cpQuery = CM.find({'is_gold' : True,'pert_id':{'$in':list(cpSet)}}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
        toDataFrame=True)
#which of the targets are lm/bing
targetList = []
for brd in cpQuery['pert_id']:
    targetList.append(compoundsExtended[brd])
cpQuery['gene_target'] = targetList

mc = mu.MongoContainer()
lmGenes = mc.gene_info.find({'is_lm':True},
        {'pr_gene_symbol':True,'is_lm':True},
        toDataFrame=True)
lmList = lmGenes['pr_gene_symbol']
cpQuery['']

# what are the LM drug-gene connections we'd expect to see?

# Which of the targets have KD?
CM = mu.CMapMongo()
goldQuery = CM.find({'is_gold' : True,'pert_desc':{'$in':geneTargets},'pert_type':'trt_sh'}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
        toDataFrame=True)
KDset = set(goldQuery['pert_iname'])
#which of the KD set are LM
mc = mu.MongoContainer()
lmStatus = mc.gene_info.find({'pr_gene_symbol':{'$in':list(KDset)}}, #, 
        {'pr_gene_symbol':True,'is_lm':True},
        toDataFrame=True)
grppedLM = lmStatus.groupby('is_lm')
iLM = grppedLM.groups[True]
kdLM = lmStatus.ix[iLM]['pr_gene_symbol']
#limit KD signatures to LM
iLM2 = goldQuery['pert_iname'].isin(kdLM)
kdLMfrm = goldQuery.ix[iLM2]

# identify which compounds target a LM with KD
dbSer = pd.Series(dbDict)
dbSer.name = 'compound_list'
cpMatch = dbSer.reindex(kdLM)
compoundsExtended = {}
for gene in cpMatch.index:
    cps = cpMatch.ix[gene]
    for cp in cps:
        compoundsExtended[cp] = gene
        #what if a compound targets multiple genes?
# get signature ids for compounds 
CM = mu.CMapMongo()
cpQuery = CM.find({'is_gold' : True,'pert_id':{'$in':compoundsExtended.keys()},},
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
        toDataFrame=True)
targetList = []
for brd in cpQuery['pert_id']:
    targetList.append(compoundsExtended[brd])
cpQuery['gene_target'] = targetList
targetGrped = cpQuery.groupby(['gene_target','cell_id'])
targetCpCounts = targetGrped.size()
targetCpCounts.name = 'target_compounds_per_cell'
tccFrame = pd.DataFrame(targetCpCounts)
#look up same info for KDs
kdGrped = kdLMfrm.groupby(['pert_iname','cell_id'])
targetKdCounts = kdGrped.size()
targetKdCounts.name = 'KDs_per_cell_line'
tkcFrame = pd.DataFrame(targetKdCounts)
#where do the cell lines overlap - KDs and Cps that target the same gene?
KdCpFrame = pd.concat([tccFrame,tkcFrame],axis=1)
KdCpFrame['KDs_per_cell_line'].notnull
nullFrm = pd.core.common.notnull(KdCpFrame)
nullVec = nullFrm['target_compounds_per_cell'] - nullFrm['KDs_per_cell_line']
KdCpFrame = KdCpFrame[~nullVec]
KdCpFrame.to_csv(wkdir+'/KD_and_compounds_targeting_the_same_gene.txt',sep='\t')


# do they have overlap of the drug in the same cell line?



### repeat for bing space

