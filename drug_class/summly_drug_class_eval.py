import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd
import matplotlib
import pandas as pd

wkdir = '/xchip/cogs/projects/cp_class_analysis/3Sept'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

# drugFile = '/xchip/cogs/projects/target_id/ctd2_annots/ctd2_merged_mapped_genes.txt'
drugFile = '/xchip/cogs/projects/cp_annot/drug_classes_AS.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')

# set up dictionary of compound classes
grpSet = set(drugLabels['CLASS'])
grpToCp = {}
for grp in grpSet:
    grpPerts = drugLabels['BRD'][drugLabels['CLASS'] == grp]
    grpToCp[grp] = list(grpPerts.values)
# compound to group dict
cpToGrp = {}
for ibrd, brd in enumerate(drugLabels['BRD']):
    cpToGrp[brd] = drugLabels['CLASS'][ibrd]

# look up brds for unkown
CM = mu.CMapMongo()
for ibrd, brd in enumerate(drugLabels['BRD']):
    if brd == '-666':
        iname = drugLabels['PERT_INAME'][ibrd]
        #look up the brd using iname
        trueBRD = CM.find({'pert_iname':iname},{'pert_id':True},limit=1)
        if trueBRD:
            trueBRD = trueBRD[0]
            #replace in dataframe
            drugLabels['BRD'][ibrd] = trueBRD

### for each drug class:
# 1) list each pert_id and iname in the group
# 2) count the number of sig sig_ids
# 3) count the number of cell lines
allBrds = drugLabels['BRD']
brdNoNull = allBrds[allBrds != '-666'] #skip -666 values
#  search for all instances
# sigFrame = CM.find({'pert_id':{'$in':list(brdNoNull)}},
#                     {'pert_id':True,'sig_id':True,'pert_iname':True,'cell_id':True},
#                     toDataFrame=True)
# is_gold and core cell line instances 
coreCells = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
sigFrame = CM.find({'pert_id':{'$in':list(brdNoNull)},'is_gold':True,'cell_id':{'$in':coreCells}},
                    {'pert_id':True,'sig_id':True,'pert_iname':True,'cell_id':True},
                    toDataFrame=True)
sigGroup = sigFrame.groupby('pert_id')
sigCounts = sigGroup.count()
sigCountSer = sigCounts['pert_id']
inameS = sigGroup.last()['pert_iname']

#get number of cell lines for each brd
nCellDict = {}
sigGrpDict = sigGroup.groups
for grp in sigGrpDict:
    iGrp = sigGrpDict[grp]
    cellSer = sigFrame.ix[iGrp]['cell_id']
    cellSet = set(cellSer)
    n_cell = len(cellSet)
    nCellDict[grp] = n_cell 
nCellSer = pd.Series(nCellDict)

classGrpBy = drugLabels.groupby('BRD')
classGrps = classGrpBy.groups
# for each compound pick one class that it belongs to
# for each compound list every group it belongs to
classDict = {}
for brd in classGrps:
    ibrd = classGrps[brd]
    grpSer = drugLabels['CLASS'][ibrd].values
    grpList = list(set(grpSer))
    classDict[brd] = grpList
sigClass = pd.Series(classDict)

### make frame of all relavant info 
# df - all drug class per cell
dfSummary = pd.concat([inameS,sigClass,sigCountSer,nCellSer], 
                keys=['iname','class','sigID_counts','n_cell'], axis=1)
# # df - one arbitrary unique class for each drug
# classFirst = classGrpBy.first()['CLASS']
# dfSummary = pd.concat([inameS,classFirst,sigCountSer,nCellSer], 
#                 keys=['iname','class','sigID_counts','n_cell'], axis=1)
sumFile = wkdir + '/class_summary.txt'
dfSummary.to_csv(sumFile)
sigFile = wkdir + '/sig_id_table.txt'
sigCellGroup = sigFrame.groupby(['pert_id','cell_id'])
sigTable = sigCellGroup.first() # pick only one sig_id for each cell line
sigTable.to_csv(sigFile)


# grpSummary = dfSummary.groupby(['class','iname']).first()
# sumFile = wkdir + '/class_summary.txt'
# grpSummary.to_csv(sumFile)

# table for rajiv's criteria
clGrpBy = drugLabels.groupby('CLASS')
clGrps = clGrpBy.groups
classDF= pd.DataFrame()
for g in clGrps:
    ig = clGrps[g]
    for ibrd in ig:
        brd = drugLabels['BRD'][ibrd]
        if brd not in sigCounts.index:
            continue
        sigsT = sigTable.ix[brd]
        sigsT['pert_id'] = brd
        sigsT['class'] = g
        sigsT['n_total_sigs'] = sigCountSer.ix[brd]
        sigsT['n_cells'] = nCellDict[brd]
        classDF = pd.concat([classDF,sigsT],axis=0)
sigFile = wkdir + '/sig_id_table.txt'
classDF.to_csv(sigFile)
