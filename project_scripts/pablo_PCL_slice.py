'''
-make an expression matrix of the top intra-connecting PCL groups
-use this for Pablo's analysis

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

wkdir = '/xchip/cogs/projects/pharm_class/svm_pcla_classifier_NOV1'
#make pso object
pso = psc.svm_pcla(out=wkdir)
self=pso

## pick 5 groups - best inter-connectors
# testGroups = ['Histone_deacetylase_1-Inhibitor',
#               'Glucocorticoid_receptor-Agonist',
#               'Proto-oncogene_tyrosine-protein_kinase_ABL1-Inhibitor',
#               'Phosphatidylinositol-4,5-bisphosphate_3-kinase_catalytic_subunit,_delta_isoform-Inhibitor',
#               '3-hydroxy-3-methylglutaryl-coenzyme_A_reductase-Inhibitor']
# load in top groups
llo = ldc.label_loader()
self.pclDict = llo.load_TTD()
#load pcl rankpoint file 
rnkpt_med_file = '/xchip/cogs/projects/pharm_class/TTd_Oct29/PCL_group_rankpt_medians.txt'
groupMedians = pd.io.parsers.read_csv(rnkpt_med_file,sep='\t')
groupMedians = groupMedians.sort('median_rankpt',ascending=False)
# make sure compounds are not counted mroe than once in a dictionary:
extendedCompoundList = []
reducedPCLDict = {}
for key in groupMedians['PCL_group']:
    value = self.pclDict[key]
    for brd in value:
        if brd in extendedCompoundList:
            value.remove(brd)
    reducedPCLDict[key] = value
    extendedCompoundList.extend(value)
self.pclDict = reducedPCLDict
n_groups = 7
testGroups = groupMedians['PCL_group'][:n_groups].values
#
brdAllGroups = []
for group in testGroups:
    brdAllGroups.extend(self.pclDict[group])
#
cellLine = 'PC3'
goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'cell_id':cellLine,'pert_dose':{'$gt':1}}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
        toDataFrame=True)
goldQuery.index = goldQuery['sig_id']
cut_by = 'pert_iname'

# asign labels
goldQuery = set_class_labels(testGroups,goldQuery,self.pclDict)

### leave only 1 or two signatures for each compound ### 
nKeep = 2
grpedBRD = goldQuery.groupby(cut_by)
keepList = []
# keep only n instances of each compound
for brd in grpedBRD.groups:
    sigs = grpedBRD.groups[brd]
    keepList.extend(sigs[:nKeep])
reducedSigFrm = goldQuery.reindex(index=keepList)
outF = wkdir + '/top_intra_connecting_compound_classes.txt'
reducedSigFrm.to_csv(outF)
# grped2 = reducedSigFrm.groupby('pert_iname')
# grped2.size()


### read in signatures ###
### write to file ####
sigList = reducedSigFrm['sig_id'].values
### load in expression data for the two sets of signatures
afPath = cmap.score_path
gt = gct.GCT()
gt.read(src=afPath,cid=sigList,rid='lm_epsilon')

outGCT = wkdir + '/top_intra_connecting_compound_classes'
gt.write(outGCT,mode='gct')
zFrm = gt.frame
# zFrm = zFrm.T
# probeIDs = zFrm.columns
# ## merge data with 
# zFrm = pd.concat([zFrm,droppedQ],axis=1)

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
