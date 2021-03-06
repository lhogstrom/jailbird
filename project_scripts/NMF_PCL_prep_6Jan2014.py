'''
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
for cellLine in cellList:
    # cellLine = 'A375'
    wkdir = '/xchip/cogs/projects/NMF/clique_n69_bing/' + cellLine
    if not os.path.exists(wkdir):
        os.mkdir(wkdir)
    # get signature annotations from cmap database
    CM = mu.CMapMongo()
    goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'cell_id':cellLine,'pert_dose':{'$gt':1}}, #, 
            {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
            toDataFrame=True)
    indRep = [x.replace(":",".") for x in goldQuery['sig_id']]
    indRep = [x.replace("-",".") for x in indRep]
    goldQuery.index = indRep
    # goldQuery.index = goldQuery['sig_id']
    dmsoQuery = CM.find({'pert_iname':'DMSO','cell_id':cellLine,'distil_nsample':{'$gte':2,'$lt':5}}, #, 
            {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
            toDataFrame=True)
    indRep = [x.replace(":",".") for x in dmsoQuery['sig_id']]
    indRep = [x.replace("-",".") for x in indRep]
    dmsoQuery.index = indRep
    # dmsoQuery['sig_id'] = indRep
    dmsoQuery['pcl_name'] = 'DMSO'
    dmsoQuery['labels'] = 99
    goldQuery = set_class_labels(testGroups,goldQuery,pclDict)
    #combine dmsos with drug perturbations
    nDMSOs = dmsoQuery.shape[0]
    iRandDmso = np.random.choice(nDMSOs-1,50,replace=False)
    goldQuery = pd.concat([goldQuery,dmsoQuery.ix[iRandDmso]],axis=0)
    ### leave only 1 or two signatures for each compound ### 
    nKeep = 2
    cut_by = 'pert_iname'
    grpedBRD = goldQuery.groupby(cut_by)
    keepList = []
    # keep only n instances of each compound
    for brd in grpedBRD.groups:
        sigs = grpedBRD.groups[brd]
        if brd == 'DMSO':
            keepList.extend(sigs) # keep all DMSO sigs
        else:
            keepList.extend(sigs[:nKeep])
    reducedSigFrm = goldQuery.reindex(index=keepList)
    outF = wkdir + '/' + cellLine + '_clique_compound_classes.v2.txt'
    reducedSigFrm.to_csv(outF,sep='\t',header=False)
    ### read in signatures ###
    ### write to file ####
    sigList = reducedSigFrm['sig_id'].values
    ### load in expression data for the two sets of signatures
    afPath = cmap.score_path
    gt = gct.GCT()
    # gt.read(src=afPath,cid=sigList,rid='lm_epsilon') # lm
    gt.read(src=afPath,cid=sigList,rid='bing') # bing
    outGCT = wkdir + '/' + cellLine + '_clique_compound_classes'
    gt.write(outGCT,mode='gctx')
    zFrm = gt.frame

# convert gctx to gct
#use java-1.7
# convert gctx to gct so it can be read by R "convert-dataset -i MCF7_top_intra_connecting_compound_classes_n130x978.gctx"
# cmd1 = 'use Java-1.7'
# os.system(cmd1)
# globRes = glob.glob(outGCT+'*.gctx')
# outF = wkdir + '/' + cellLine + '_top_intra_connecting_compound_classes.v2.txt'
# outF = cellLine + '_top_intra_connecting_compound_classes.v2.txt'
# cmd2 = 'convert-dataset -i ' + globRes[0]
# os.system(cmd2)

# ### load in Pablo's dir
# cellLine = 'A375'
# wkdir = '/xchip/cogs/projects/NMF/clique_n69/' + cellLine
# Hfile = '/xchip/cogs/projects/NMF/clique_n69/A375/A375_clique_compound_classes_n328x978.H.k9.txt'
# WFile = '/xchip/cogs/projects/NMF/clique_n69/A375/A375_clique_compound_classes_n328x978.W.k9.txt'
# aFile = '/xchip/cogs/projects/NMF/clique_n69/A375/A375_clique_compound_classes.v2.txt'
cellLine = 'A549'
wkdir = '/xchip/cogs/projects/NMF/clique_n69/' + cellLine
Hfile = '/xchip/cogs/projects/NMF/clique_n69/A549/A549_clique_compound_classes_n436x978.H.k9.txt'
WFile = '/xchip/cogs/projects/NMF/clique_n69/A549/A549_clique_compound_classes_n436x978.W.k9.txt'
aFile = '/xchip/cogs/projects/NMF/clique_n69/A549/A549_clique_compound_classes.v2.txt'
Hmtrx = pd.io.parsers.read_csv(Hfile,sep='\t',index_col=0) #,header=True
Hmtrx = Hmtrx.T
Wmtrx = pd.io.parsers.read_csv(WFile,sep='\t',index_col=0) #,header=True
anntFrm = pd.read_csv(aFile,sep='\t',header=False,index_col=0)
anntFrm.columns = reducedSigFrm.columns
groupSer = pd.Series(index=anntFrm.index,data=anntFrm.pcl_name)
if (groupSer.index == Hmtrx.index).all():
    iZip = zip(*[groupSer.values,groupSer.index.values])
    mCol = pd.MultiIndex.from_tuples(iZip, names=['pcl_name','sig_id'])
    Hmtrx.index = mCol
pclGrped = Hmtrx.groupby(level='pcl_name')
graphDir = wkdir + '/graphs'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
for grp in pclGrped.groups.keys():
    grpH = pclGrped.get_group(grp)
    plt.imshow(grpH,
        interpolation='nearest',
        cmap=cm.gray_r)
    ytcks = list(grpH.index.get_level_values('sig_id'))
    xtcks = list(grpH.columns)
    plt.xticks(np.arange(len(xtcks)), xtcks,rotation=75)
    plt.yticks(np.arange(len(ytcks)),ytcks)
    plt.colorbar()
    grpMod = grpMod = ''.join(e for e in grp if e.isalnum())
    outF = os.path.join(graphDir,grpMod+'.png')
    plt.savefig(outF, bbox_inches='tight')
    plt.close()


# #assume sig by component
# ### run classifier on the 
# predictDict = {}
# for sig in Hmtrx.index:
#     droppedFrm = Hmtrx[Hmtrx.index != sig] # remove test signature from training
#     trainFrm = droppedFrm.reindex(columns=probeIDs)
#     labelsTrain = droppedFrm['labels'].values
#     C = 1.0  # SVM regularization parameter
#     svc = svm.SVC(kernel='linear', C=C).fit(trainFrm.values, labelsTrain)
#     zTest = zFrm.ix[sig,probeIDs]
#     linPred = svc.predict(zTest.values)

# ### run sig_introspect on the same signatures
# outIntrospect = '/xchip/cogs/hogstrom/analysis/pablos_NMF_analysis/' + cellLine + 'sig_introspect_by_cell'
# if not os.path.exists(outIntrospect):
#     os.mkdir(outIntrospect)
# qSer = reducedSigFrm['sig_id']
# os.chdir(wkdir)
# qSer.to_csv(outF,index=False,header=False)
# #run sig_introspect
# cmd = ' '.join(['rum -q hour',
#      '-d sulfur_io=100',
#      '-o ' + outIntrospect,
#      '-x sig_introspect_tool ',
#      '--sig_id ' + outF,
#      '--query_group cell_id',
#      '--metric wtcs',
#      '--out ' + outIntrospect])
# os.system(cmd)

# #load in sig_introspect result
# file_rnkpt = '/xchip/cogs/hogstrom/analysis/pablos_NMF_analysis/MCF7sig_introspect_by_cell/nov12/my_analysis.sig_introspect_tool.4658128.0/self_rankpt_n79x79.gctx'
# gt = gct.GCT(file_rnkpt)
# gt.read()
# rnkptFrm = gt.frame

# groupSigs = rnkptFrm.index.values
# #get info for these signatures
# CM = mu.CMapMongo()
# goldQuery = CM.find({'is_gold' : True,'sig_id':{'$in':list(groupSigs)}}, #
#         {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
#         toDataFrame=True)
# goldQuery.index = goldQuery['sig_id']
# goldQuery = set_class_labels(testGroups,goldQuery,self.pclDict)

# # make individual heatmaps for each compound class - recorded in an
# # individual cell line
# cut_by = 'pcl_name'
# grpedBRD = goldQuery.groupby(cut_by)
# for group in grpedBRD.groups:
#     sigs = grpedBRD.groups[group]
#     tmpRnkpt = rnkptFrm.reindex(index=sigs,columns=sigs)
#     ### heatmap code - by group
#     fig = plt.figure(1, figsize=(10, 8))
#     plt.suptitle(group + ' compound group',fontsize=14, fontweight='bold')
#     plt.subplot(111)
#     plt.title('mean_rankpt_4')
#     colors.set_color_map()
#     plt.imshow(tmpRnkpt.values,
#             interpolation='nearest',
#             vmin=-100, 
#             vmax=100)
#     ytcks = [goldQuery['pert_iname'][x] for x in sigs]
#     plt.xticks(np.arange(len(ytcks)), ytcks,rotation=75)
#     plt.yticks(np.arange(len(ytcks)),ytcks)
#     plt.colorbar()
#     out = outIntrospect + '/' + group + '_rnkpt_heatmap.png'
#     fig.savefig(out, bbox_inches='tight')
#     plt.close()



# ### psudeo code for baseline expression heatmaps:
# # has dave done this before with pert_explorer?
# # Just extending this to more drug-gene relationships?
# # how to run sig_introspect on things that are 

# # log2 abundance (baseline expression) --> y-axis is differential
# # gene expression

# # screen for instances were drug-gene relationships are higher in 
# # are higher in cell lines with baseline expression of a given gene

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

