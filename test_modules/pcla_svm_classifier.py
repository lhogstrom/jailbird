'''
PCLA (pharmacological class analyzer) classifier

Run svm classification on various PCLA drug classes

Larson Hogstrom, 9/2013
'''
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
import cmap.util.progress as update
import multiprocessing
import time

class svm_pcla(object):
    '''
    Class to serve PCLA classification
    Parameters
    ----------
    out : str
        out directory path
    '''
    def __init__(self,
                out):
        '''
        Initialize a new instance of classifier
        
        '''
        # set output directories
        self.out = out
        if not os.path.exists(self.out):
            os.mkdir(self.out)
        # set core cell lines
        coreCells = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
        self.core_cell_lines = coreCells

    def set_classes(self):
        '''
        specify source of class labels
        Parameters
        ----------
        '''
        ### load in data for individual groups
        llo = ldc.label_loader()
        self.pclDict = llo.load_TTD()
        ## pick 5 groups - best inter-connectors
        testGroups = ['Histone_deacetylase_1-Inhibitor',
                      'Glucocorticoid_receptor-Agonist',
                      'Proto-oncogene_tyrosine-protein_kinase_ABL1-Inhibitor',
                      'Phosphatidylinositol-4,5-bisphosphate_3-kinase_catalytic_subunit,_delta_isoform-Inhibitor',
                      '3-hydroxy-3-methylglutaryl-coenzyme_A_reductase-Inhibitor']
        brdAllGroups = []
        for group in testGroups:
            brdAllGroups.extend(self.pclDict[group])
        self.all_group_cps = brdAllGroups
        self.test_groups = testGroups

    def set_clique_n69(self):
        '''
        load 69 PCLs currated by Rajiv 

        Parameters
        ----------
        '''
        ### load in data for individual groups
        llo = ldc.label_loader()
        self.pclDict = llo.load_clique_set_n69()
        self.all_group_cps = self.pclDict.values()
        self.test_groups = self.pclDict.keys()

    def PCL_vs_DMSO(self,max_signatures_per_cp=20,n_test_max=False):
        '''
        -grab equal amounts of DMSO and signatures from a PCL class
        -test one PCL at a time

        Parameters
        ----------
        n_test_max : int
            -max number of PCL groups to incorporate into the classifier 
            -if set to False, all groups are tested
        '''
        for group_name in self.test_groups:
            group_cps = self.pclDict[group_name]
            CM = mu.CMapMongo()
            # set minimum dose
            cpQuery = CM.find({'is_gold' : True,'pert_id':{'$in':group_cps},'pert_dose':{'$gt':1}}, #, 
                    {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
                    toDataFrame=True)
            # inameGrped = cpQuery.groupby('pert_iname')
            cpQuery.index = cpQuery['sig_id']
            cpQuery = self.set_class_labels(cpQuery)
            droppedQ = self.cut_signatures(cpQuery,nKeep=max_signatures_per_cp,cut_by='pert_iname')
            droppedGrped = droppedQ.groupby('pert_iname')
            droppedGrped.size()
        # if groups_to_model == None:
        #     groups_to_model = self.pclDict.keys()
        # brdAllGroups = []
        # for group in groups_to_model:
        #     brdAllGroups.extend(self.pclDict[group])
        # CM = mu.CMapMongo()
        # # set minimum dose
        # goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'pert_dose':{'$gt':1}}, #, 
        #         {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
        #         toDataFrame=True)
        # goldQuery.index = goldQuery['sig_id']
        # # asign drug class labels
        # goldQuery = self.set_class_labels(goldQuery)
        # # reduce signatures to prevent overfitting to one compound
        # droppedQ = self.cut_signatures(goldQuery,nKeep=max_signatures_per_cp)
        # sigList = droppedQ['sig_id'].values
        # ### load in expression data for the two sets of signatures
        # afPath = cmap.score_path
        # gt = gct.GCT()
        # gt.read(src=afPath,cid=sigList,rid='lm_epsilon')
        # zFrm = gt.frame
        # zFrm = zFrm.T
        # probeIDs = zFrm.columns
        # ## merge data with 
        # zFrm = pd.concat([zFrm,droppedQ],axis=1)

    def test_classes_incrementally(self,rnkpt_med_file,n_test_max=False):
        '''
        -start from the most internally consistent PCL - and move down the list 
        -incrementally increase the number of groups added to the classifier

        Parameters
        ----------
        rnkpt_med_file : str
            path to a file containing the median summly rankpoint values for each
            group (output from the pcla tool)
        n_test_max : int
            -max number of PCL groups to incorporate into the classifier 
            -if set to False, all groups are tested
        '''
        ### load in data for individual groups
        llo = ldc.label_loader()
        self.pclDict = llo.load_TTD()
        #load pcl rankpoint file 
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
        # set incrament of groups
        if n_test_max:
            max_groups = n_test_max
        else:
            max_groups = groupMedians.shape[0]
        group_range = np.arange(2,max_groups+1)
        n_group_accuracy = {}
        for n_groups in group_range:
            print "testing " + str(n_groups) + " number of classes"
            testGroups = groupMedians['PCL_group'][:n_groups].values
            self.test_groups = testGroups
            self.classification_across_cell(groups_to_model=testGroups,loo_type='by_cp',max_signatures_per_cp=3)
            n_group_accuracy[n_groups] = self.model_accuracy_across_cells
        self.n_group_accuracy = n_group_accuracy

    def incremental_group_size_graph():
        '''
        -Graph the results for building SVMs of multiple numbers of groups
        -XY scatter

        Parameters
        ----------
        none

        ''' 
        groupSer = pd.Series(self.n_group_accuracy)
        plt.plot(groupSer.index.values, groupSer.values, 'o-')
        plt.ylim([0,1])
        plt.xlabel('PCLs classified')
        plt.ylabel('classification accuracy')
        plt.title('top most intra-connected PCL groups')
        outF = os.path.join(self.out,'vary_PCL_groups_classfied.png')
        plt.savefig(outF, bbox_inches='tight')
        plt.close()

    def classification_by_cell(self,loo_type='by_cp'):
        '''
        -For each of the specified cell lines, build a separate classifier
        -evaluate model with leave one out cross val.
        
        Parameters
        ----------
        loo_type : str
            strategy for leave one out validation:
                'by_cp' - leaves out all signatures for a given compounds
                'by_sig' - leaves out individual signatures 
        '''        
        combinedFrm = pd.DataFrame()
        accuracyDict = {}
        for cellLine in self.core_cell_lines:
            CM = mu.CMapMongo()
            # goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'cell_id':cellLine}, #, 
            #         {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
            #         toDataFrame=True)
            # set minimum dose
            goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':self.all_group_cps},'cell_id':cellLine,'pert_dose':{'$gt':1}}, #, 
                    {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
                    toDataFrame=True)
            goldQuery.index = goldQuery['sig_id']
            # asign drug class labels
            goldQuery = self.set_class_labels(goldQuery)
            # reduce signatures to prevent overfitting to one compound
            droppedQ = self.cut_signatures(goldQuery)
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
            if loo_type == 'by_cp':
                zFrm['svm_prediction'] = np.nan
                cpSet = set(zFrm['pert_id'])
                # loop through the compounds - leave out in building the model then test
                for brd in cpSet:
                    brd_match = zFrm['pert_id'] == brd
                    droppedFrm = zFrm[~brd_match] # remove test signature from training
                    trainFrm = droppedFrm.reindex(columns=probeIDs)
                    labelsTrain = droppedFrm['labels'].values
                    C = 1.0  # SVM regularization parameter
                    svc = svm.SVC(kernel='linear', C=C).fit(trainFrm.values, labelsTrain)
                    zTest = zFrm.ix[brd_match,probeIDs]
                    linPred = svc.predict(zTest.values)
                    zFrm['svm_prediction'][zTest.index] = linPred
            if loo_type == 'by_sig':
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
            self.modelFrame = combinedFrm
            self.model_accuracy = accuracyDict

    def classification_across_cell(self,loo_type='by_cp',max_signatures_per_cp=3,groups_to_model=None):
        '''
        -build a single classifier treating observations from different
        cell lines equally
        -evaluate model with leave one out cross val.
        
        Parameters
        ----------
        groups_to_model : list
            -list of group names in the pclDict
            -default is to use all keys
        loo_type : str
            strategy for leave one out validation:
                'by_cp' - leaves out all signatures for a given compounds
                'by_sig' - leaves out individual signatures 
        max_signatures_per_cp : int
            maximum number of signatures per compound to incorporate into the classifier
            (to avoid overfitting to compounds with many signatures)
        '''        
        #start a update indicator
        progress_bar = update.DeterminateProgressBar('SVM calculation')
        if groups_to_model == None:
            groups_to_model = self.pclDict.keys()
        brdAllGroups = []
        for group in groups_to_model:
            brdAllGroups.extend(self.pclDict[group])
        CM = mu.CMapMongo()
        # set minimum dose
        goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'pert_dose':{'$gt':1}}, #, 
                {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
                toDataFrame=True)
        goldQuery.index = goldQuery['sig_id']
        # asign drug class labels
        goldQuery = self.set_class_labels(goldQuery)
        # reduce signatures to prevent overfitting to one compound
        droppedQ = self.cut_signatures(goldQuery,nKeep=max_signatures_per_cp)
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
        if loo_type == 'by_cp':
            zFrm['svm_prediction'] = np.nan
            cpSet = set(zFrm['pert_id'])
            tupList = [(zFrm,brd,probeIDs) for brd in cpSet]
            # run SVM in parallel
            prog = update.DeterminateProgressBar('self-connection graph builder')
            n_procs=8
            pool = multiprocessing.Pool(n_procs)
            rs = pool.map_async(_svm_worker,tupList)
            pool.close() # No more work
            while (True):
                if (rs.ready()): break
                remaining = rs._number_left
                prog.show_message('SVM evaluation - {0} tasks to complete'.format(remaining))
                time.sleep(0.1)
            results = rs.get()
            predictedSer = pd.Series()
            for result in results:
                predictedSer = predictedSer.append(result)
            zFrm['svm_prediction'] = predictedSer
        if loo_type == 'by_sig':
            predictDict = {}
            for ii,sig in enumerate(zFrm.index):
                progress_bar.update('running SVM and signature validation - ' + sig, ii, len(zFrm.index))
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
        accuracyArray = zFrm['labels'] == zFrm['svm_prediction']
        accuracyRate = accuracyArray.sum()/float(accuracyArray.shape[0])
        self.model_accuracy_across_cells = accuracyRate
        self.signature_frame = zFrm

    def set_class_labels(self,sigInfoFrm):
        '''
        set known labels for test data

        Parameters
        ----------
        sigInfoFrm : pandas dataFrame
            dataFrame of signature info where index are sig_ids
        ''' 
        sigInfoFrm['labels'] = np.nan
        sigInfoFrm['pcl_name'] = 'null'
        for igroup,group in enumerate(self.test_groups):
            grpMembers = self.pclDict[group]
            iMatch = sigInfoFrm['pert_id'].isin(grpMembers)
            sigInfoFrm['labels'][iMatch] = igroup
            sigInfoFrm['pcl_name'][iMatch] = group
        return sigInfoFrm

    def cut_signatures(self,sigInfoFrm,nKeep=2,cut_by='pert_id'):
        '''
        limit the number signatures to prevent over fitting to a single compound

        Parameters
        ----------
        sigInfoFrm : pandas dataFrame
            dataFrame of signature info where index are sig_ids
        nKeep : int
            number of signatures to keep for each compound
        cut_by : str
            sig_info field to group and cut by
        
        Returns
        ----------
        reducedSigFrm : pandas dataFrame
            sigInfoFrm with less signatures - about even for each compound
            dataFrame of signature info where index are sig_ids

        ''' 
        grpedBRD = sigInfoFrm.groupby(cut_by)
        keepList = []
        # keep only n instances of each compound
        for brd in grpedBRD.groups:
            sigs = grpedBRD.groups[brd]
            keepList.extend(sigs[:nKeep])
        reducedSigFrm = sigInfoFrm.reindex(index=keepList)
        # grped = reducedSigFrm.groupby('pcl_name')
        # grped.size()
        return reducedSigFrm

def _svm_worker(argTup):
    '''
    worker to build SVM model and validation for one drug

    Parameters
    ----------
    argTup : tuple
        contains this set of arguments:
            zFrm : pandas dataFrame
                dataFrame of signature info z score expression data
            brd : str
                string of a pert_id
            probeIDs : list or index of probeIDs
                probe space to perform calculations (eg, all landmark genes)
    Returns
    ----------
    predSer : pandas Series
        index = sig_ids
        values = predicted drug class
    ''' 
    zFrm = argTup[0]
    brd = argTup[1]
    probeIDs = argTup[2]
    brd_match = zFrm['pert_id'] == brd
    droppedFrm = zFrm[~brd_match] # remove test signature from training
    trainFrm = droppedFrm.reindex(columns=probeIDs)
    labelsTrain = droppedFrm['labels'].values
    C = 1.0  # SVM regularization parameter
    svc = svm.SVC(kernel='linear', C=C).fit(trainFrm.values, labelsTrain)
    zTest = zFrm.ix[brd_match,probeIDs]
    linPred = svc.predict(zTest.values)
    predSer = pd.Series(linPred,index=zTest.index)
    predSer.index.name = brd
    return predSer


