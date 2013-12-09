'''
-load in NMF results
-perform SVM on H matrix
-compare to SVM built on raw data

Larson Hogstrom, 11/2013
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
import pandas as pd
import cmap
import os
import cmap.analytics.pcla as pcla
import cmap.plot.colors as colors
from sklearn import svm, datasets


### load in Pablo's dir
resDir = '/xchip/cogs/projects/NMF/lincs_core_cell_lines/A375'
W9file = 'A375_top_intra_connecting_compound_classes_n109x978.W.k9.gct'
H9file = 'A375_top_intra_connecting_compound_classes_n109x978.H.k9.gct'
annotFile = 'A375_top_intra_connecting_compound_classes.v2.txt'
# w9 = gct.GCT('/'.join([resDir,W9file]))
# w9.read()
annotFrm = pd.io.parsers.read_csv('/'.join([resDir,annotFile]),sep='\t',index_col=0,header=None)

Hmtrx = pd.io.parsers.read_csv('/'.join([resDir,H9file]),sep='\t',skiprows=1,header=1,index_col=0) #,header=True, index_col=0
Hmtrx = Hmtrx.drop('Description',1)
sigs = Hmtrx.columns.values
Hmtrx = Hmtrx.T

annotFrm = annotFrm.reindex(list(sigs))
labels = annotFrm.ix[:,5]
sigIDs = list(annotFrm.ix[:,10].values)

# run SVM 
# train with whole dataset
# svc = svm.SVC(kernel='linear', C=C).fit(Hmtrx.values, labels)
# linPred = svc.predict(Hmtrx.values)
#leave one out
predictDict = {}
for sig in Hmtrx.index:
    trainFrm = Hmtrx[Hmtrx.index != sig] # remove test signature from training
    labelsTrain = labels[labels.index != sig]
    # trainFrm = droppedFrm.reindex(columns=probeIDs)
    C = 1.0  # SVM regularization parameter
    svc = svm.SVC(kernel='linear', C=C,probability=True).fit(trainFrm.values, labelsTrain.values)
    zTest = Hmtrx.ix[sig,:]
    linPred = svc.predict(zTest.values)
    linPred_prob = svc.predict_proba(zTest.values)
    linPred_df = svc.decision_function(zTest.values)
    predictDict[sig] = linPred[0]
predSer = pd.Series(predictDict)
predSer.name = 'svm_prediction'
predSer = predSer.reindex(labels.index)
accuracyArray = predSer.values == labels.values
accuracyRate = accuracyArray.sum()/float(accuracyArray.shape[0])

# get expression values for sigs
afPath = cmap.score_path
gt = gct.GCT()
gt.read(src=afPath,cid=sigIDs,rid='lm_epsilon')
zFrm = gt.frame
zFrm = zFrm.True

#compare SVM with raw data
#replace index
if (zFrm.index == annotFrm.ix[labels.index,10]).all():
    zFrm.index = annotFrm.index
predictDict = {}
for sig in zFrm.index:
    trainFrm = zFrm[zFrm.index != sig] # remove test signature from training
    labelsTrain = labels[labels.index != sig]
    C = 1.0  # SVM regularization parameter
    svc = svm.SVC(kernel='linear', C=C).fit(trainFrm.values, labelsTrain.values)
    zTest = zFrm.ix[sig,:]
    linPred = svc.predict(zTest.values)
    predictDict[sig] = linPred[0]
predSer = pd.Series(predictDict)
predSer.name = 'svm_prediction'
predSer = predSer.reindex(labels.index)
accuracyArray = predSer.values == labels.values
accuracyRate = accuracyArray.sum()/float(accuracyArray.shape[0])


# componentIndex = Hmtrx.columns
# ## add drug labels
# labelFrm = reducedSigFrm.ix[:,['labels','pcl_name']]
# # labelFrm = labelFrm.reindex(Hmtrx.index)
# Hmtrx = pd.concat([Hmtrx,labelFrm],axis=0,ignore_index=False)

