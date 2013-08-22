#! /usr/bin/env python
'make a roc curve to look at drug-gene connnections from summly results'

import numpy as np
import pylab as pl
from sklearn import svm, datasets
from sklearn.utils import shuffle
from sklearn.metrics import roc_curve, auc

import os
import cmap.util.mongo_utils as mu
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib

### load in CTD2 summly results 
work_dir = '/xchip/cogs/hogstrom/analysis/summly/drug_gene_roc'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

drugFile = '/xchip/cogs/projects/target_id/ctd2_annots/ctd2_merged_mapped_genes.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
targetSet = set(drugLabels['gene_dw_annot'])
### pair brd and pert_iname
allCtd2 = drugLabels['pert_id']
allinames = drugLabels['pert_iname']
pDescDict = {}
for ibrd,brd in enumerate(allCtd2):
    pDescDict[brd] = allinames[ibrd]

####################################
### set up compound annotations ####
####################################

### make target_dict
targetSheetF = '/xchip/cogs/projects/target_id/7June2014/A2_DrugBank_targets_tab.txt'
targetDict = {}
pDescDict = {}
with open(targetSheetF,'rt') as f:
    for string in f:
        splt = string.split('\r')
        for i,line in enumerate(splt):
            splt2 = line.split('\t')
            pID = splt2[0] #the pert_id listed the line
            pDesc = splt2[1]
            targets = splt2[2:]
            targets = [x for x in targets if x != '']
            # targets = targets.split(';')
            if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
                continue
            else:
                targetDict[pID] = targets
                pDescDict[pID] = pDesc
targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_drug_targets.txt'
# targetDict = {}
# pDescDict = {}
with open(targetSheetF,'rt') as f:
    for string in f:
        splt = string.split('\r')
        for i,line in enumerate(splt):
            splt2 = line.split('\t')
            pID = splt2[0] #the pert_id listed the line
            pDesc = splt2[1]
            targets = splt2[2]
            targets = targets.split(';')
            targets = [x for x in targets if x != '']
            if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
                continue
            else:
                targetDict[pID] = targets
                pDescDict[pID] = pDesc

########################################
###### full ctd2 summly results ######
########################################

brdGrpList = []
grpSet = set(drugLabels['gene_dw_annot'])
# group to compound dict
grpToCp = {}
for grp in grpSet:
    grpPerts = drugLabels['pert_id'][drugLabels['gene_dw_annot'] == grp]
    grpToCp[grp] = list(grpPerts.values)
    brdGrpList.extend(grpPerts.values)
# compound to group dict
cpToGrp = {}
for ibrd, brd in enumerate(drugLabels['pert_id']):
    cpToGrp[brd] = drugLabels['gene_dw_annot'][ibrd]
# get list of cps in summly dir
# basePath = work_dir + '/ctd2_sig_query'
basePath = '/xchip/cogs/hogstrom/analysis/summly/cp_class/ctd2_sig_query'
dateID = 'aug01/my_analysis.sig_summly_tool.2013080119394091'
# dateID = 'jul29/my_analysis.sig_summly_tool.2013072914520591'
summDir = '/'.join([basePath,dateID])
cpDirs = [f for f in os.listdir(summDir) if os.path.isdir(summDir+'/'+f)]

sumScoreFrm = pd.DataFrame()
for ibrd,brd in enumerate(cpDirs):    
    inFile = '/'.join([summDir,
                    brd,
                    brd+'_summly.txt'])
    sumRes = pd.io.parsers.read_csv(inFile,sep='\t')
    # filter to only cps / cgs
    pd.io.parsers.read_csv
    cpRes = sumRes[sumRes['pert_type'] == 'trt_cp']
    cpRes['rank'] = np.arange(1,len(cpRes)+1)
    cgsRes = sumRes[sumRes['pert_type'] == 'trt_sh.cgs']
    cgsRes['rank'] = np.arange(1,len(cgsRes)+1)
    oeRes = sumRes[sumRes['pert_type'] == 'trt_oe']
    oeRes['rank'] = np.arange(1,len(oeRes)+1)
    #sort by summly results
    tmpFrame = pd.DataFrame(data=oeRes['sum_score'].values,index=oeRes['pert_iname'],columns=[brd])
    sumScoreFrm = pd.concat([sumScoreFrm,tmpFrame],axis=1)

### create a matrix of expected relationships
dgMtrx = sumScoreFrm.copy()
dgMtrx.ix[:,:] = np.zeros_like(sumScoreFrm)
for brd in dgMtrx.columns:
    if not brd in targetDict:
        print brd + ' not in target dictionary'
        continue
    targets = targetDict[brd]
    if targets:
        for target in targets:
            if not target in dgMtrx.index:
                print target + ' not in target dictionary'
                continue
            dgMtrx.ix[target,brd] = 1

### ROC analysis for target-ID classes ###
sumScores = sumScoreFrm.unstack()
dg_labels = pd.DataFrame(dgMtrx.unstack(),columns=['label'])
# labelFrm = pd.DataFrame(data=sumScores.values,index=sumScores.index,columns=[brd])
labelFrm = pd.DataFrame(data=sumScores,columns=['sum_scores'])
labelFrm = pd.concat([labelFrm,dg_labels],axis=1)
labelSrt = labelFrm[sumScores.notnull()]
labelSrt = labelSrt.sort(columns='sum_scores' ,ascending=False)

# fpr, tpr, thresholds = roc_curve(dg_labels, sumScores)
fpr, tpr, thresholds = roc_curve(labelSrt['label'].values, labelSrt['sum_scores'].values)
roc_auc = auc(fpr, tpr)
print "Area under the ROC curve : %f" % roc_auc

# Plot ROC curve
pl.clf()
pl.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
pl.plot([0, 1], [0, 1], 'k--')
pl.xlim([0.0, 1.0])
pl.ylim([0.0, 1.0])
pl.xlabel('False Positive Rate')
pl.ylabel('True Positive Rate')
pl.title('Expected Drug-OE ROC')
pl.legend(loc="lower right")
pl.show()

### ROC analysis for functional classes ###
sumScoreFrm = pd.DataFrame()
for ibrd,brd in enumerate(cpDirs):    
    inFile = '/'.join([summDir,
                    brd,
                    brd+'_summly.txt'])
    sumRes = pd.io.parsers.read_csv(inFile,sep='\t')
    # filter to only cps / cgs
    pd.io.parsers.read_csv
    cpRes = sumRes[sumRes['pert_type'] == 'trt_cp']
    cpRes['rank'] = np.arange(1,len(cpRes)+1)
    cgsRes = sumRes[sumRes['pert_type'] == 'trt_sh.cgs']
    cgsRes['rank'] = np.arange(1,len(cgsRes)+1)
    oeRes = sumRes[sumRes['pert_type'] == 'trt_oe']
    oeRes['rank'] = np.arange(1,len(oeRes)+1)
    #sort by summly results
    tmpFrame = pd.DataFrame(data=cpRes['sum_score'].values,index=cpRes['pert_id'],columns=[brd])
    sumScoreFrm = pd.concat([sumScoreFrm,tmpFrame],axis=1)

### create a matrix of expected relationships
grpMtrx = sumScoreFrm.copy()
grpMtrx.ix[:,:] = np.zeros_like(sumScoreFrm)
for brd in grpMtrx.columns:
    grp = cpToGrp[brd]
    if grp == '-666':
        continue
    grpBrds = grpToCp[grp]
    grpBrds.remove(brd) #take out self from group list
    for brd_conn in grpBrds:
        if not brd_conn in grpMtrx.index:
            print brd_conn + ' not in target dictionary'
            continue
        grpMtrx.ix[brd_conn,brd] = 1

### 
sumScores = sumScoreFrm.unstack()
dg_labels = pd.DataFrame(grpMtrx.unstack(),columns=['label'])
labelFrm = pd.DataFrame(data=sumScores,columns=['sum_scores'])
labelFrm = pd.concat([labelFrm,dg_labels],axis=1)
labelSrt = labelFrm[sumScores.notnull()]
labelSrt = labelSrt.sort(columns='sum_scores' ,ascending=False)
# roc calculation
fpr, tpr, thresholds = roc_curve(labelSrt['label'].values, labelSrt['sum_scores'].values)
roc_auc = auc(fpr, tpr)
print "Area under the ROC curve : %f" % roc_auc
# Plot ROC curve
pl.clf()
pl.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
pl.plot([0, 1], [0, 1], 'k--')
pl.xlim([0.0, 1.0])
pl.ylim([0.0, 1.0])
pl.xlabel('False Positive Rate')
pl.ylabel('True Positive Rate')
pl.title('expected drug-drug ROC')
pl.legend(loc="lower right")
pl.show()




### examine each functional gorup
for grpGene in grpToCp:
    if grpGene == '-666':
        continue
    grp = grpToCp[grpGene]
    if not grp: # skip if grp is empty
        continue
    # matrices for grp connections
    nGrp = len(grp)
    grp_sum_score = np.zeros((nGrp,nGrp))
    grp_PercSummly = np.zeros((nGrp,nGrp))
    grp_rank = np.zeros((nGrp,nGrp))
    for ibrd,brd in enumerate(grp):


    if not brd in targetDict:
        print brd + ' not in target dictionary'
        continue
    targets = targetDict[brd]
    if targets:
        for target in targets:
            if not target in grpMtrx.index:
                print target + ' not in target dictionary'
                continue
            grpMtrx.ix[target,brd] = 1


### scratch 
# tmpFrame = pd.DataFrame(data=cgsRes['sum_score'].values,index=cgsRes['pert_iname'],columns=[brd])
# inameSer = cgsRes['pert_iname'].copy()
# inameSer.sort()
# #series
# tmpSer = pd.Series(data=cgsRes['sum_score'].values,index=cgsRes['pert_iname'],name=brd)
# tempSer = tmpSer.reindex(index=inameSer)
# sumScoreFrm = sumScoreFrm.append(tmpSer)

# sumScoreFrm = sumScoreFrm.concat(tmpSer)
# sumScoreFrm = pd.concat([sumScoreFrm,tmpFrame],axis=1)
