'''
-prep CMAP matrices for NMF analysis

Larson Hogstrom, 04/2014
'''
import numpy as np
import pylab as pl
from matplotlib import cm
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import cmap.io.gmt as gmt
import pandas as pd
import cmap
import os
import cmap.analytics.pcla as pcla
import cmap.plot.colors as colors
import glob
import shutil

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

cellList = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines

### set parameters
probeSpace = 'lm_epsilon' # lm_epsilon or bing
nDMSO = 50
nKeep = 2 # number of signatures per drug
for cell in cellList:
    print(cell)
    prefix = cell + '_c20_' + probeSpace
    wkdir = '/xchip/cogs/projects/NMF/cpd_groups_n147/' + prefix
    if not os.path.exists(wkdir):
        os.mkdir(wkdir)
    #load in clique annotations and matrix
    # cFile = '/xchip/cogs/sig_tools/sig_cliquescore_tool/sample/cp_clique_n69/clique.gmt'
    # cFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140221/cliques.gmt'
    cFile = '/xchip/cogs/projects/pharm_class/rnwork/cliques/cpd_groups_n147.gmt'
    cliqueGMT = gmt.read(cFile)
    cliqFrm = pd.DataFrame(cliqueGMT)
    # set grouping structures 
    pclDict = {}
    for x in cliqFrm.iterrows():
        pclDict[x[1]['id']] = set(x[1]['sig'])
    # create list of all compounds
    brdAllGroups = []
    for group in pclDict:
        brdAllGroups.extend(pclDict[group])
    brdAllGroups.append('DMSO')
    brdAllGroups = list(set(brdAllGroups))
    testGroups = cliqFrm['id'].values
    # extract signatures and expression data for every group member
    ### Write a big matrix for all cell lines ####
    # get signature annotations from cmap database
    CM = mu.CMapMongo()
    goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'cell_id':cell,'pert_dose':{'$gt':1}}, #, 
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
        # reducedSigFrm = pd.concat([reducedSigFrm,goldCell],axis=0)
        #combine dmsos with drug perturbations
        dmsoCell = cellDMSOGrped.get_group(cell)
        nDMSOs = dmsoCell.shape[0]
        iRandDmso = np.random.choice(nDMSOs-1,nDMSO,replace=False)
        reducedSigFrm = pd.concat([reducedSigFrm,goldCell,dmsoQuery.ix[iRandDmso]],axis=0) 
    reducedSigFrm.index.name = 'mod_sig_id'
    outF = wkdir + '/clique_compound_classes.v2.txt'
    # reducedSigFrm.to_csv(outF,sep='\t',header=False) #without header
    reducedSigFrm.to_csv(outF,sep='\t',header=True) #with header
    ### read in signatures ###
    sigList = list(set(reducedSigFrm['sig_id'].values))
    ### load in expression data for the two sets of signatures
    afPath = cmap.score_path
    gt = gct.GCT()
    gt.read(src=afPath,cid=sigList,rid=probeSpace) # lm
    # gt.read(src=afPath,cid=sigList,rid='bing') # bing
    ### write to file ####
    outGCT = wkdir + '/clique_compound_classes'
    gt.write(outGCT,mode='gctx')

### convert GCTX files to GCT ### 
#use java-1.7
probeSpace = 'lm_epsilon' # lm_epsilon or bing
for cell in cellList:
    print(cell)
    prefix = cell + '_c20_' + probeSpace
    wkdir = '/xchip/cogs/projects/NMF/cpd_groups_n147/' + prefix
    os.chdir(wkdir)
    # cmd1 = 'use Java-1.7'
    # os.system(cmd1)
    outGCT = wkdir + '/clique_compound_classes'
    globRes = glob.glob(outGCT+'*.gctx')
    print(globRes[0])
    cmd2 = 'convert-dataset -i ' + globRes[0]
    os.system(cmd2)
#     dir1 = glob.glob(wkdir+'/*02')
#     dir2 = glob.glob(dir1[0]+'/*')
#     dir3 = glob.glob(dir2[0]+'/*')
#     gctFile = dir3[0]
#     shutil.copy(gctFile, wkdir)
