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

wkdir = '/xchip/cogs/projects/NMF/clique_n69_all_cell_lines'
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
for cell in cellGoldGrped.groups.keys():
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

### Analyze W matrix ###
wkdir = '/xchip/cogs/projects/NMF/clique_n69_all_cell_lines'
Hfile = '/xchip/cogs/projects/NMF/clique_n69_all_cell_lines/clique_compound_classes_n2873x978.H.k9.txt'
WFile = '/xchip/cogs/projects/NMF/clique_n69_all_cell_lines/clique_compound_classes_n2873x978.W.k9.gct'
aFile = '/xchip/cogs/projects/NMF/clique_n69_all_cell_lines/clique_compound_classes.v2.txt'
Hmtrx = pd.io.parsers.read_csv(Hfile,sep='\t',index_col=0) #,header=True
Hmtrx = Hmtrx.T
Wmtrx = pd.io.parsers.read_csv(WFile,sep='\t',index_col=0) #,header=True
anntFrm = pd.read_csv(aFile,sep='\t',header=False,index_col=0)
anntFrm.columns = reducedSigFrm.columns
# groupSer = pd.Series(index=anntFrm.index,data=anntFrm.pcl_name)
# if (groupSer.index == Hmtrx.index).all():
#     iZip = zip(*[groupSer.values,groupSer.index.values])
#     mCol = pd.MultiIndex.from_tuples(iZip, names=['pcl_name','sig_id'])
#     Hmtrx.index = mCol
# pclGrped = Hmtrx.groupby(level='pcl_name')
graphDir = wkdir + '/graphs_uniform_max_sort'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
# flH = np.float64(Hmtrx.values)
maxVal = Hmtrx.max(axis=1).max()
for r in cliqFrm.iterrows():
    grp = r[1]['id']
    brds = r[1]['sig']
    anntMtch = anntFrm[anntFrm.pert_id.isin(brds)]
    grpH = Hmtrx.reindex(anntMtch.index)
    meanVec = grpH.describe().ix['mean']
    ### correlate each observed sigature with mean vector
    # corrMtrx = np.corrcoef(meanVec,grpH)
    # corrVec = corrMtrx[:-1,-1]
    # iSort = np.argsort(corrVec)
    # grpH = grpH.ix[iSort,:] # sort acording to corr with mean
    ### take top three components - order acording to their strenght
    iTop3 = meanVec.order(ascending=False).index[:3]
    sortedTop = grpH.ix[:,iTop3].sort()
    topSum = sortedTop.sum(axis=1).order(ascending=False)
    grpH = grpH.ix[topSum.index,:] # sort acording to corr with mean
    Hfloat = np.float64(grpH.values)
    fig = plt.figure(figsize=(40, 20), dpi=50)
    plt.imshow(Hfloat,
        interpolation='nearest',
        cmap=cm.gray_r,
        vmax=maxVal)
    ytcks = list(grpH.index)
    xtcks = list(grpH.columns)
    plt.xticks(np.arange(len(xtcks)), xtcks,rotation=75)
    plt.yticks(np.arange(len(ytcks)),ytcks)
    plt.colorbar()
    grpMod = grpMod = ''.join(e for e in grp if e.isalnum())
    outF = os.path.join(graphDir,grpMod+'.png')
    plt.savefig(outF, bbox_inches='tight')
    plt.close()


