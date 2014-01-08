'''
-Make one large data matrix for all cell lines
-make an expression matrix of the top intra-connecting PCL groups
-use this for Pablo's NMF analysis

Larson Hogstrom, 01/2014
'''
import test_modules.pcla_svm_classifier as psc
import numpy as np
import pylab as pl
from sklearn import svm, datasets
from matplotlib import cm
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import test_modules.load_TTD_drug_class as ldc
import cmap.io.gct as gct
import cmap.io.gmt as gmt
import pandas as pd
import cmap
import os
import cmap.analytics.pcla as pcla
import cmap.plot.colors as colors
import glob

def set_class_labels(test_groups,sigInfoFrm,pclDict):
    '''
    set known labels for test data

    Parameters
    ----------
    sigInfoFrm : pandas dataFrame
        dataFrame of signature info where index are sig_ids
    ''' 
    sigInfoFrm['labels'] = np.nan
    sigInfoFrm['pcl_name'] = 'null'
    for igroup,group in enumerate(test_groups):
        grpMembers = pclDict[group]
        iMatch = sigInfoFrm['pert_id'].isin(grpMembers)
        sigInfoFrm['labels'][iMatch] = igroup
        sigInfoFrm['pcl_name'][iMatch] = group
    return sigInfoFrm

wkdir = '/xchip/cogs/projects/NMF/clique_n69_bing'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

#load in clique annotations and matrix
cFile = '/xchip/cogs/sig_tools/sig_cliquescore_tool/sample/cp_clique_n69/clique.gmt'
cliqueGMT = gmt.read(cFile)
cliqFrm = pd.DataFrame(cliqueGMT)

# set grouping structures 
pclDict = {}
for x in cliqFrm.iterrows():
    pclDict[x[1]['id']] = set(x[1]['sig'])
#
brdAllGroups = []
for group in pclDict:
    brdAllGroups.extend(pclDict[group])
brdAllGroups.append('DMSO')
brdAllGroups = list(set(brdAllGroups))
testGroups = cliqFrm['id'].values

# extract signatures and expression data for every group member
# for each cell line
cellList = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
### Write a big matrix for all cell lines ####
wkdir = '/xchip/cogs/projects/NMF/clique_n69_all_cell_lines'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
# get signature annotations from cmap database
CM = mu.CMapMongo()
goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'cell_id':{'$in':cellList},'pert_dose':{'$gt':1}}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
        toDataFrame=True)
indRep = [x.replace(":",".") for x in goldQuery['sig_id']]
indRep = [x.replace("-",".") for x in indRep]
goldQuery.index = indRep
# goldQuery.index = goldQuery['sig_id']
dmsoQuery = CM.find({'pert_iname':'DMSO','cell_id':{'$in':cellList},'distil_nsample':{'$gte':2,'$lt':5}}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
        toDataFrame=True)
indRep = [x.replace(":",".") for x in dmsoQuery['sig_id']]
indRep = [x.replace("-",".") for x in indRep]
dmsoQuery.index = indRep
# dmsoQuery['sig_id'] = indRep
dmsoQuery['pcl_name'] = 'DMSO'
dmsoQuery['labels'] = 99
goldQuery = set_class_labels(testGroups,goldQuery,pclDict)
# group by cell line
cellGoldGrped = goldQuery.groupby('cell_id')
cellDMSOGrped = dmsoQuery.groupby('cell_id')
reducedSigFrm = pd.DataFrame()
for cell in cellGrped.groups.keys():
    goldCell = cellGoldGrped.get_group(cell)
    ### leave only 1 or two signatures for each compound ### 
    nKeep = 2
    cut_by = 'pert_iname'
    grpedBRD = goldCell.groupby(cut_by)
    keepList = []
    # keep only n instances of each compound
    for brd in grpedBRD.groups:
        sigs = grpedBRD.groups[brd]
        if brd == 'DMSO':
            keepList.extend(sigs) # keep all DMSO sigs
        else:
            keepList.extend(sigs[:nKeep])
    goldCell = goldCell.reindex(index=keepList)
    # compounds only
    reducedSigFrm = pd.concat([reducedSigFrm,goldCell],axis=0)
    #combine dmsos with drug perturbations
    # dmsoCell = cellDMSOGrped.get_group(cell)
    # nDMSOs = dmsoCell.shape[0]
    # iRandDmso = np.random.choice(nDMSOs-1,50,replace=False)
    # reducedSigFrm = pd.concat([reducedSigFrm,goldCell,dmsoQuery.ix[iRandDmso]],axis=0)
outF = wkdir + '/clique_compound_classes.v2.txt'
reducedSigFrm.to_csv(outF,sep='\t',header=False)
### read in signatures ###
### write to file ####
sigList = list(set(reducedSigFrm['sig_id'].values))
### load in expression data for the two sets of signatures
afPath = cmap.score_path
gt = gct.GCT()
gt.read(src=afPath,cid=sigList,rid='lm_epsilon') # lm
# gt.read(src=afPath,cid=sigList,rid='bing') # bing
outGCT = wkdir + '/clique_compound_classes'
gt.write(outGCT,mode='gctx')
zFrm = gt.frame

