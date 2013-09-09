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

wkdir = '/xchip/cogs/sig_tools/sig_summly/pcl'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

# drugFile = '/xchip/cogs/projects/target_id/ctd2_annots/ctd2_merged_mapped_genes.txt'
# drugFile = '/xchip/cogs/projects/cp_annot/drug_classes_AS.txt'
drugFile = '/xchip/cogs/sig_tools/sig_summly/pcl/pcl_classes.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
#repalce ugly characters
drugLabels['class'].str.replace("/","-")
drugLabels['class'].str.replace(" ","_")
drugLabels['class'].str.replace("/","_")
drugLabels['class'].str.replace("&","_")
drugLabels['class'].str.replace("?","_")
drugLabels['class'].str.replace("(","_")
drugLabels['class'].str.replace(")","_")

# set up dictionary of compound classes
grpSet = set(drugLabels['class'])
grpToCp = {}
for grp in grpSet:
    grpPerts = drugLabels['pert_id'][drugLabels['class'] == grp]
    grpToCp[grp] = list(grpPerts.values)
# compound to group dict
cpToGrp = {}
for ibrd, brd in enumerate(drugLabels['pert_id']):
    cpToGrp[brd] = drugLabels['class'][ibrd]

# look up brds for unkown
CM = mu.CMapMongo()
for ibrd, brd in enumerate(drugLabels['pert_id']):
    if brd == '-666':
        iname = drugLabels['PERT_INAME'][ibrd]
        #look up the brd using iname
        trueBRD = CM.find({'pert_iname':iname},{'pert_id':True},limit=1)
        if trueBRD:
            trueBRD = trueBRD[0]
            #replace in dataframe
            drugLabels['pert_id'][ibrd] = trueBRD

### for each drug class:
# 1) list each pert_id and iname in the group
# 2) count the number of sig sig_ids
# 3) count the number of cell lines
allBrds = drugLabels['pert_id']
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

classGrpBy = drugLabels.groupby('pert_id')
classGrps = classGrpBy.groups
# for each compound pick one class that it belongs to
# for each compound list every group it belongs to
classDict = {}
for brd in classGrps:
    ibrd = classGrps[brd]
    grpSer = drugLabels['class'][ibrd].values
    grpList = list(set(grpSer))
    classDict[brd] = grpList
sigClass = pd.Series(classDict)

### make frame of all relavant info 
# df - all drug class per cell
dfSummary = pd.concat([inameS,sigClass,sigCountSer,nCellSer], 
                keys=['iname','class','sigID_counts','n_cell'], axis=1)
# # df - one arbitrary unique class for each drug
# classFirst = classGrpBy.first()['class']
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
clGrpBy = drugLabels.groupby('class')
clGrps = clGrpBy.groups
classDF= pd.DataFrame()
for g in clGrps:
    ig = clGrps[g]
    for ibrd in ig:
        brd = drugLabels['pert_id'][ibrd]
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

# load in Rajiv's file
pclFile = '/xchip/cogs/sig_tools/sig_summly/pcl/pcl_sig_info.txt'
pclTbl = pd.read_csv(pclFile,sep='\t')
sumSpaceFile = '/xchip/cogs/sig_tools/sig_query/results/a2_internal_wtcs.lm/query_info_n73597.txt'
summSpace = pd.read_csv(sumSpaceFile,sep='\t')

#sort based on pert_id
pclGrp = pclTbl.groupby('pert_id')
pclDict= pclGrp.groups
for brd in pclDict:
    tmpTbl = pclTbl.ix[pclDict[brd]] #table for the one pert_id
    tmpSer = tmpTbl['sig_id'] # series of sig_ids
    brdFile = wkdir + '/list/' + brd + '.grp'
    tmpSer.to_csv(brdFile,index=False) # write sigs for each brd to a file
#make dictionary of inames
firstTbl = pclGrp.first()
inameDict = {}
for brd in firstTbl.index:
    inameDict[brd] = firstTbl.ix[brd]['pert_iname']

#submit summly jobs to lsf - CELL LINE MATCH MODE
for brd in pclDict:
    summlyMtrx = '/xchip/cogs/sig_tools/sig_query/results/a2_internal_wtcs.lm/'
    outDir = '/xchip/cogs/sig_tools/sig_summly/pcl/summly_out_match'
    querySpace = wkdir + '/' + brd + '.grp'
    cmd = ' '.join(['rum -q hour -x sig_summly_tool',
             summlyMtrx,
             '--query_space ' + querySpace,
             '--group_query true',
             '--out ' + outDir])
    os.system(cmd)

#submit summly jobs to lsf - CELL LINE INDEPENDENT MODE
for brd in pclDict:
    summlyMtrx = '/xchip/cogs/sig_tools/sig_query/results/a2_internal_wtcs.lm/'
    outDir = '/xchip/cogs/sig_tools/sig_summly/pcl/summly_out_independent'
    querySpace = wkdir + '/' + brd + '.grp'
    cmd = ' '.join(['rum -q hour -x sig_summly_tool',
             summlyMtrx,
             '--query_space ' + querySpace,
             '--group_query false',
             '--out ' + outDir])
    os.system(cmd)

########################################
###### full ctd2 summly results ######
########################################

# take the average of a pairwise matrix
def av_mtrx(mtrx):
    nm = len(mtrx)
    avMtrx = np.zeros((nm,nm))
    for i1 in range(nm):
        for i2 in range(nm):
            val1 = mtrx[i1,i2]
            val2 = mtrx[i2,i1]
            avMtrx[i1,i2] = np.mean([val1,val2])
            # avMtrxUp = np.triu(avMtrx,k=1)
            iUp = np.tril_indices(nm)
            avMtrx[iUp] = np.nan
    return avMtrx

# get list of cps in summly dir
# basePath = work_dir + '/sig_query'
basePath = '/xchip/cogs/sig_tools/sig_summly/pcl/summly_out/sep06'
summDir = [f for f in os.listdir(basePath) if os.path.isdir(basePath+'/'+f)]
summPath = [basePath+'/'+f for f in summDir]
#put the path for each cp in a dictionary
cpPathDict = {}
for p in summPath:
    cpDirs = [f for f in os.listdir(p) if os.path.isdir(p+'/'+f)]
    cp = cpDirs[0]
    cpPathDict[cp] = p + '/' + cp
# check to make sure all cps have summly result
summlyLeftOut = []
for brd in pclDict:
    if not brd in cpPathDict:
        summlyLeftOut.appen(brd)

graphDir = wkdir + '/graphs_Sept9'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)

### examine each functional gorup
sumScoreDict = {} #matrix of connections for each group
percSummlyDict = {}
for grpName in grpToCp:
    if grpName == '-666':
        continue
    grp = grpToCp[grpName]
    grp = [cp for cp in grp if cp in cpPathDict] # leave out compounds that don't have summly data
    if not grp: # skip if grp is empty
        continue
    # matrices for grp connections
    nGrp = len(grp)
    grp_sum_score = np.zeros((nGrp,nGrp))
    grp_PercSummly = np.zeros((nGrp,nGrp))
    grp_rank = np.zeros((nGrp,nGrp))
    for ibrd,brd in enumerate(grp):
        if not brd in cpPathDict:
            continue
        inFile = '/'.join([cpPathDict[brd],
                        brd+'_summly.txt'])
        sumRes = pd.io.parsers.read_csv(inFile,sep='\t')
        # filter to only cps / cgs
        cpRes = sumRes[sumRes['pert_type'] == 'trt_cp']
        cpRes['rank'] = np.arange(1,len(cpRes)+1)
        cgsRes = sumRes[sumRes['pert_type'] == 'trt_sh.cgs']
        cgsRes['rank'] = np.arange(1,len(cgsRes)+1)
        oeRes = sumRes[sumRes['pert_type'] == 'trt_oe']
        oeRes['rank'] = np.arange(1,len(oeRes)+1)
        #check group connection 
        for ibrd2, brd2 in enumerate(grp):
            indSum = cpRes[cpRes['pert_id'] == brd2]['sum_score_4']
            if not indSum:
                print brd + ' ' + brd2 + ' not compared' 
                continue
            sumScore = indSum.values[0]
            indrank = cpRes[cpRes['pert_id'] == brd2]['rank']
            rank = indrank.values[0]
            percSummly = rank / float(len(cpRes))
            grp_sum_score[ibrd,ibrd2] = sumScore
            grp_PercSummly[ibrd,ibrd2] = percSummly
            grp_rank[ibrd,ibrd2] = rank
    #make pairwise matrix into a pandas dataframe
    sumScoreFrm = pd.DataFrame(grp_sum_score,index=grp,columns=grp)


    #take averages
    av_grp_sum_score = av_mtrx(grp_sum_score)
    av_grp_PercSummly = av_mtrx(grp_PercSummly)
    av_grp_rank = av_mtrx(grp_rank)
    # store averages
    sumScoreDict[grpName] = grp_sum_score
    percSummlyDict[grpName] = grp_PercSummly
    ### print group heatmap
    fig = plt.figure(1, figsize=(20, 8))
    plt.suptitle(grpName + ' compound group',fontsize=14, fontweight='bold')
    plt.subplot(121)
    plt.title('percent summly rank')
    plt.imshow(av_grp_PercSummly,
            interpolation='nearest',
            cmap=matplotlib.cm.Greens_r,
            vmin=0, 
            vmax=1)
    ytcks = [inameDict[x] for x in grp]
    plt.xticks(np.arange(len(grp)), ytcks,rotation=75)
    plt.yticks(np.arange(len(grp)),ytcks)
    plt.colorbar()
    plt.subplot(122)
    plt.title('sum_score')
    plt.imshow(av_grp_sum_score,
            interpolation='nearest',
            cmap=matplotlib.cm.RdBu_r,
            vmin=-1, 
            vmax=1)
    plt.xticks(np.arange(len(grp)), ytcks,rotation=75)
    plt.yticks(np.arange(len(grp)),ytcks)
    plt.colorbar()
    outF = os.path.join(graphDir,grpName + '_compound_group_heatmap.png')
    fig.savefig(outF, bbox_inches='tight')
    plt.close()

#make boxplot of all connections
sumScoreList = []
percSumList = []
tickList = []
for gName in sumScoreDict:
    # sum score setup
    m1 = sumScoreDict[gName]
    upMtrx = av_mtrx(m1)
    flatM = upMtrx.flatten()
    flatM = flatM[~np.isnan(flatM)] # remove nan
    sumScoreList.append(flatM)
    # percent summly setup
    m2 = percSummlyDict[gName]
    upMtrx2 = av_mtrx(m2)
    flatM2 = upMtrx2.flatten()
    flatM2 = flatM2[~np.isnan(flatM2)] # remove nan
    percSumList.append(flatM2)
    #names
    tickList.append(gName)
plt.boxplot(sumScoreList)
plt.xticks(np.arange(len(tickList)),tickList,rotation=45)
plt.xlabel('compound class - by Molecular target',fontweight='bold')
plt.ylabel('sum_score',fontweight='bold')
plt.title('distribution of sum scores for CTD2 compound class',fontweight='bold')

plt.boxplot(percSumList)
plt.xticks(np.arange(len(tickList)),tickList,rotation=45)
plt.xlabel('compound class - by molecular target',fontweight='bold')
plt.ylabel('percent summly',fontweight='bold')
plt.title('distribution of summly percents for CTD2 compound class',fontweight='bold')


