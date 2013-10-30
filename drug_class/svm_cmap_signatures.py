import numpy as np
import pylab as pl
from sklearn import svm, datasets
from matplotlib import cm
import cmap.util.mongo_utils as mu

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
realVPred = zFrm.ix[:,['labels','svm_prediction']]
accuracyArray = zFrm['labels'] == zFrm['svm_prediction']
accuracyRate = accuracyArray.sum()/float(accuracyArray.shape[0])


### first iteration with two common drugs
CM = mu.CMapMongo()
goldQuery = CM.find({'is_gold' : True}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
        toDataFrame=True)
grped = goldQuery.groupby(['pert_id','cell_id'])
grpSize = grped.size()
grpSize = grpSize.order(ascending=False)
grp1 = grpSize.index[0]
grp2 = grpSize.index[1]
igrp1 = grped.groups[grp1]
igrp2 = grped.groups[grp2]
grp1Quer = goldQuery.ix[igrp1]
grp2Quer = goldQuery.ix[igrp2]
grp1Sigs = grp1Quer['sig_id']
grp2Sigs = grp2Quer['sig_id']
sigList = list(grp1Sigs)
sigList.extend(grp2Sigs)
# load in expression data 
afPath = cmap.score_path
gt = gct.GCT()
gt.read(src=afPath,cid=sigList,rid='lm_epsilon')
zFrm = gt.frame
zFrm = zFrm.T
probeIDs = zFrm.columns
## asign labels
pert_ids = [x.split(':')[1][:13] for x in zFrm.index]
# create labels
labels = []
for brd in pert_ids:
  if brd == 'BRD-A19500257':
    labels.append(1)
  else:
    labels.append(0)
labels = np.array(labels)
zFrm['labels'] = labels