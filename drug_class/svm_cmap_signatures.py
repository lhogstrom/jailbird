import numpy as np
import pylab as pl
from sklearn import svm, datasets
from matplotlib import cm
import cmap.util.mongo_utils as mu
import test_modules.load_TTD_drug_class as ldc
import cmap.io.gct as gct
import pandas as pd
import cmap

### current parameters:
# 1) cell lines
# 2) test groups
# 3) number of sigs per compound
# 4) is_gold vs. non (dose)
# 5) svm model type: linear, poly, etc
# 6) C = 1.0  # SVM regularization parameter
# 7) compounds that belong to more than one group - how to 
# deal with these?

# Things to try:
# 1) Testing model on extra signatures (non LOA) - exta signatures from mongo
# 2) Testing classifier on DOS compounds 
# 3) build classifiers for each PCL individually (group vs. DMSO) - which PCL works best?

### load in data for individual groups
llo = ldc.label_loader()
pclDict = llo.load_TTD()
## pick 5 groups - best inter-connectors
testGroups = ['Histone_deacetylase_1-Inhibitor',
              'Glucocorticoid_receptor-Agonist',
              'Proto-oncogene_tyrosine-protein_kinase_ABL1-Inhibitor',
              'Phosphatidylinositol-4,5-bisphosphate_3-kinase_catalytic_subunit,_delta_isoform-Inhibitor',
              '3-hydroxy-3-methylglutaryl-coenzyme_A_reductase-Inhibitor']
brdAllGroups = []
for group in testGroups:
  brdAllGroups.extend(pclDict[group])
### look up data for the compounds
coreCells = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
combinedFrm = pd.DataFrame()
accuracyDict = {}
for cellLine in coreCells:
    CM = mu.CMapMongo()
    # goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'cell_id':cellLine}, #, 
    #         {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
    #         toDataFrame=True)
    # set minimum dose
    goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'cell_id':cellLine,'pert_dose':{'$gt':1}}, #, 
            {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
            toDataFrame=True)
    goldQuery.index = goldQuery['sig_id']
    #make labels
    goldQuery['labels'] = np.nan
    goldQuery['pcl_name'] = 'null'
    for igroup,group in enumerate(testGroups):
        grpMembers = pclDict[group]
        iMatch = goldQuery['pert_id'].isin(grpMembers)
        goldQuery['labels'][iMatch] = igroup
        goldQuery['pcl_name'][iMatch] = group
    # count group members
    grpedBRD = goldQuery.groupby('pert_id')
    keepList = []
    # keep only n instances of each compound
    nKeep = 2
    for brd in grpedBRD.groups:
        sigs = grpedBRD.groups[brd]
        keepList.extend(sigs[:nKeep])
    droppedQ = goldQuery.reindex(index=keepList)
    grped = droppedQ.groupby('pcl_name')
    grped.size()
    sigList = droppedQ['sig_id'].values
    ### load in expression data for the two sets of signatures
    afPath = cmap.score_path
    gt = gct.GCT()
    gt.read(src=afPath,cid=sigList,rid='lm_epsilon')
    zFrm = gt.frame
    zFrm = zFrm.T
    probeIDs = zFrm.columns
    ## merge data with 
    zFrm = pd.concat([zFrm,droppedQ],axis=1)
    ### perform leave one out validation
    predictDict = {}
    for sig in zFrm.index:
      droppedFrm = zFrm[zFrm.index != sig] # remove test signature from training
      trainFrm = droppedFrm.reindex(columns=probeIDs)
      labelsTrain = droppedFrm['labels'].values
      C = 1.0  # SVM regularization parameter
      svc = svm.SVC(kernel='linear', C=C).fit(trainFrm.values, labelsTrain)
      zTest = zFrm.ix[sig,probeIDs]
      linPred = svc.predict(zTest.values)
      predictDict[sig] = linPred[0]
    predSer = pd.Series(predictDict)
    predSer.name = 'svm_prediction'
    zFrm = pd.concat([zFrm,pd.DataFrame(predSer)],axis=1)
    combinedFrm = pd.concat([combinedFrm,zFrm],axis=0)
    accuracyArray = zFrm['labels'] == zFrm['svm_prediction']
    accuracyRate = accuracyArray.sum()/float(accuracyArray.shape[0])
    accuracyDict[cellLine] = accuracyRate

## do LOA for all signatures of one compound at a time
## move down the list of PCLs

# split data into two and train on the rest
nSigs = combinedFrm.shape[0]
randVec = np.random.randint(2, size=nSigs)
testBool = np.ma.make_mask(randVec)
trainBool = ~testBool
trainFrm = combinedFrm[trainBool]
labelsTrain = trainFrm['labels'].values
trainFrm = trainFrm.reindex(columns=probeIDs)
testFrm = combinedFrm[testBool]
testLM = testFrm.reindex(columns=probeIDs)
C = 1.0  # SVM regularization parameter
svc = svm.SVC(kernel='linear', C=C).fit(trainFrm.values, labelsTrain)
linPred = svc.predict(testLM.values)
testFrm['svm_prediction'] = linPred
accuracyArray = testFrm['labels'] == testFrm['svm_prediction']


## accuracy across all cell lines
accuracyArray = combinedFrm['labels'] == combinedFrm['svm_prediction']
accuracyRate = accuracyArray.sum()/float(accuracyArray.shape[0])
### types of miss-classifications
wrongFrm = zFrm[~accuracyArray]
# wrongFrm['labels'] + ' ' +  wrongFrm['svm_prediction']

### descriptions of 
grped = combinedFrm.groupby(['pcl_name','cell_id'])
# grped = combinedFrm.groupby('pcl_name')
grpCounts = grped.size()

#########
## old ##
#########
## load in cmap data for svm
#pick pick a select group of compounds
drugList = ['BRD-K81418486','BRD-A19500257','BRD-A19037878','BRD-A75409952','BRD-A79768653']
cellLine = 'MCF7'
CM = mu.CMapMongo()
goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':drugList},'cell_id':cellLine}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
        toDataFrame=True)
goldQuery.index = goldQuery['sig_id']
#make labels
goldQuery['labels'] = np.nan
for ibrd,brd in enumerate(drugList):
    iMatch = goldQuery['pert_id'] == brd
    goldQuery['labels'][iMatch] = ibrd
# grped = goldQuery.groupby('pert_id')
sigList = goldQuery['sig_id'].values

#load in expression data for the two sets of signatures
afPath = cmap.score_path
gt = gct.GCT()
gt.read(src=afPath,cid=sigList,rid='lm_epsilon')
zFrm = gt.frame
zFrm = zFrm.T
probeIDs = zFrm.columns
## merge data with 
zFrm = pd.concat([zFrm,goldQuery],axis=1)

