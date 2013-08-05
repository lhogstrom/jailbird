import glob, HTML
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection as fdr
import cmap.util.mongo_utils as mu
from cmap.tools import sig_slice_tool
from cmap.io import gct,plategrp,rnk
import cmap.util.progress as progress
import subprocess
import datetime
import cmap.util.tool_ops as to
import cmap.analytics.dgo as dgo
# import cmap.analytics.dgo_oracle as dgo_oracle

work_dir = '/xchip/cogs/projects/target_id/HOG_8July2013_not_gold'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

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
# pert_List = CM.find({'sig_id':{'$regex':'DOSBIO'},'pert_iname':{'$regex':'BRD'}},{'sig_id':True,'pert_id':True,'pert_iname':True})
CM = mu.CMapMongo()
hogList = CM.find({'sig_id':{'$regex':'HOG'}},{'pert_id':True})
hogSet = set(hogList)
#check how many doses were recorded for each cp
doseDict = {}
for brd in hogSet:
    idoses = CM.find({'sig_id':{'$regex':'HOG'},'pert_id':brd},{'pert_dose':True})
    doseDict[brd] = set(idoses)
hogSet.remove('CMAP-HSF-HOGA1')
hogSet.remove('DMSO')

### retrieve targets for HOG cps 
hogTargetDict = {}
for brd in hogSet:
    if brd in targetDict:
        hogTargetDict[brd] = targetDict[brd]

####################################
### make drug-gene object ##########
####################################

dg = dgo.QueryTargetAnalysis(out=work_dir)
dg.add_dictionary(targetDict=hogTargetDict)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
dg.make_result_frames(gp_type='KD',metric='spearman')
dg.test_known_connections(gp_type='KD',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='apriori_graphs',
                        conn_thresh=.05,
                        make_graphs=True,
                        dose_graph=True,
                        n_rand=10000,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='KD',
                metric='spearman',
                outName='apriori_FDR_pass',
                alpha=0.2,
                make_graphs=True,
                specificity_graph=True)
dg.fdr_html_summary(fdrDir='test_FDR2',specificity_graph=True)    

####################################
### dgo scratch 
####################################

gp_type='KD'
metric='spearman'
pDescDict=pDescDict
outName='test_dose_graphs'
conn_thresh=.05
make_graphs=True
dose_graph=True
n_rand=10000
connection_test='two_sided'

graphDir = dg.outputdir + '/' + outName
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
#get brds from result dataframe
brdSkipped = []
cgsSkipped = []
BRDsTested = []
for ind in dg.dfRank.index:
    brd = ind[0]
    BRDsTested.append(brd)
brdRsltSet = set(BRDsTested)
#get cgs tested 
if gp_type == 'KD':
    cols = dg.dfRank.columns
if gp_type == 'OE':
    cols = [col[0] for col in dg.dfRank.columns]
countDict = {}
cellCountDict = {}
nConnectionDict = {}
connCellDict = {}
pDict = {}
pVec = []
prog = progress.DeterminateProgressBar('Connection test')
for ibrd,brd in enumerate(dg.targetDict):
    # skip pert if not in result file
    prog.update(brd,ibrd,len(dg.targetDict))
    targets = dg.targetDict[brd]
    if not brd in brdRsltSet:
        brdSkipped.append(brd)
        for target in targets:
            countDict[brd + ':' + target] = 0
        continue
    cpRes = dg.dfCS.ix[brd]
    cpRank = dg.dfRank.ix[brd]
    meanSer = cpRes.mean()
    meanRnk = cpRank.mean()
    nullCnt = pd.isnull(cpRes)
    #how many cell lines were both the pert and target tested in
    valCounts = nullCnt.shape[0] - nullCnt.sum(axis=0)
    for target in targets:
        tarList = [inst for inst in cols if inst.split('_')[0] == target]
        if len(tarList) == 0: #skip if drug target not tested
            cgsSkipped.append(target)
            continue
        for ind in tarList:
            rnkSer = cpRank[ind]
            if gp_type == 'OE':
                rnkSer = rnkSer.unstack()
            rnkSer = rnkSer[rnkSer.notnull()]
            csSer = cpRes[ind]
            if gp_type == 'OE':
                csSer = csSer.unstack()
            csSer = csSer[csSer.notnull()]
            #skip if cgs not tested in the same cell line as cp
            if len(rnkSer) == 0:
                cgsSkipped.append(ind)
                continue
            ### calculate p-value - based on percent rank products
            rnkSmll = rnkSer/100 #return rank percent betweeen 0 and 1
            testStat = rnkSmll.prod()
            n_obs = rnkSer.shape[0]
            ### simulate random draws from percent rank list
            permMtrx = np.random.rand(n_obs,n_rand)
            nullDist = permMtrx.prod(axis=0)
            #number of null values more extreme than observed 
            PosExVals = nullDist[nullDist<testStat] # positive connections
            NegExVals = nullDist[nullDist>testStat] # negative connections
            if connection_test == 'one_sided':
                nExtreme = len(PosExVals)
                pVal = (nExtreme+1)/float(len(nullDist))
            if connection_test == 'two_sided':
                NnExtreme = len(NegExVals)
                PnExtreme = len(PosExVals)
                #devide p-value by 2 since you are testing twice as many hypothesis
                if PnExtreme <= NnExtreme:
                    pVal = (PnExtreme+1)/(float(len(nullDist))/2)
                else:
                    pVal = (NnExtreme+1)/(float(len(nullDist))/2)
            pVec.append(pVal)
            pDict[brd + ':' + ind] = pVal
            countDict[brd + ':' + ind] = len(rnkSer) # number of instances for the drug-target pair
            if gp_type == 'KD':
                cells =[in1.split('_')[1] for in1 in rnkSer.index]
            if gp_type == 'OE':
                cells =[in1[0].split('_')[1] for in1 in rnkSer.index]
            cellSet = set(cells)
            cellCountDict[brd + ':' + ind] = len(cellSet)
            #describe connections which pass arbitrary threshold
            n_connections = len(rnkSer[rnkSer < (conn_thresh*100)])
            rnkSerThresh = rnkSer[rnkSer < (conn_thresh*100)]
            if rnkSerThresh.any():
                cells = [ind1.split('_')[1] for ind1 in rnkSerThresh.index]
                cells = list(set(cells))
                connCellDict[brd + ':' + ind] = (':').join(cells)
            nConnectionDict[brd + ':' + ind] = n_connections
            #make summary output
            outF = os.path.join(graphDir,brd +'_' + ind + '_drug-target_summary.txt')
            dg.make_CS_summary(brd,pDescDict[brd],rnkSer,csSer,pVal,outF,gp_type,metric)
            #best connection table
            outBestTbl = os.path.join(graphDir,'best_drug-target_summary.txt')
            dg.make_CS_summary(brd,pDescDict[brd],rnkSer,csSer,pVal,outBestTbl,gp_type,metric,best_connection_tbl=False)
            if make_graphs:
                ### cs wadden gram
                csGraph = os.path.join(graphDir,brd +'_' + ind + '_connections.png')
                dg.connection_dot_plot(csSer,
                                        pDescDict[brd],
                                        ind,
                                        csGraph,
                                        gp_type=gp_type,
                                        axis=metric)
                ### percent rank wadden gram
                csGraph = os.path.join(graphDir,brd +'_' + ind + '_percent_rank.png')
                dg.connection_dot_plot(rnkSer,
                                        pDescDict[brd],
                                        ind,
                                        csGraph,
                                        gp_type=gp_type,
                                        axis='perc_rank')
                if dose_graph:
                    dg.__connection_dose_plot(rnkSer,
                                            pDescDict[brd],
                                            ind,
                                            graphDir,
                                            gp_type=gp_type,
                                            axis='perc_rank')
                    dg.__connection_dose_plot(csSer,
                                            pDescDict[brd],
                                            ind,
                                            graphDir,
                                            gp_type=gp_type,
                                            axis=metric)  

# connection dose module

    # def __connection_dose_plot(self,inSer,pert_iname,gene1,outdir,gp_type='KD',axis='perc_rank'):
    #     '''
    #     make a dot plot (WaddenGram) for a set of dose data
    #     inputs:
    #     inSer - pandas Series of percent rank or wtcs scores between a drug-gene pair
    #     pert_iname - description of your starting drug signature
    #     gene1 - genetic perturbation you are connecting to

    #     axis:
    #     -percent rank 
    #     -wtcs
        # '''

# axis=metric
# inSer = csSer
axis='perc_rank'
inSer = rnkSer
pert_iname = pDescDict[brd]
gene1 = ind
outdir = graphDir
gp_type=gp_type

platePres = [x.split(':')[0] for x in inSer.index]
kd_tps = [x.split('_')[1]+ '_' + x.split('_')[2] for x in platePres]
doses = [float(x.split(':')[-1]) for x in inSer.index]
doseFrame = pd.DataFrame({'rank':inSer,'dose': doses,'cell_tp':kd_tps})
tpSet = set(kd_tps)
for cell_tp in tpSet:
    smFrame = doseFrame[doseFrame['cell_tp'] == cell_tp]
    smFrame = smFrame.sort('dose')
    brd = smFrame.index[0].split(':')[1][:13]
    doseRange = range(1,len(smFrame)+1)
    doseStr = [str(x) for x in smFrame['dose']]
    fig1 = plt.figure(figsize = [12,10])
    # h1 = plt.scatter(doseSpace[:18],rnks.values[:18],label='A549',color='purple',s=60,alpha=.4)
    h1 = plt.scatter(doseRange,smFrame['rank'],label=cell_tp,color='green',s=60,alpha=.4)
    plt.legend()
    if axis == 'perc_rank':
        plt.ylim((-5, 100))
    else:
        plt.ylim((-1, 1))
    plt.xlabel('dose um',fontsize=20,fontweight='bold')
    plt.xticks(range(1,10), doseStr, rotation = 45)
    plt.title(pert_iname + ' - ' + cell_tp + ' connection - ' + gp_type,fontweight='bold')
    if axis == 'perc_rank':
        plt.ylabel('percent rank',fontsize=20,fontweight='bold')
        plt.savefig(os.path.join(outdir,'_'.join([brd,cell_tp,gp_type,gene1,'dose_percent_rank.png'])))
    else:
        plt.ylabel(axis,fontsize=20,fontweight='bold')
        plt.savefig(os.path.join(outdir,'_'.join([brd,cell_tp,gp_type,gene1,'dose_metric.png'])))
    plt.close()


### 1+2 dose-consistent connection
# requires strong dose connection at 1 conecntration
# consistant directionality at two adjacent concentrations

#make sure rnkFrame is sorted by dose

#was the strongest connection positive or negative?
max1 = max(smFrame['rank'])
min1 = min(smFrame['rank'])

minThresh1 = 1
minThresh2 = 2
maxThresh1 = 85
maxThresh2 = 75

#is min val less than thresh1?
if min1 <= minThresh1:
    minFlag = 1
    iSig = smFrame['rank'].argmin() #index of most extreme
    #find indices of extreme
    iSigs = smFrame[smFrame['rank'] < minThresh1].index
    # connDirect = 'pos'
else:
    minFlag = 0

#is max val greater than thresh1?
if max1 >= maxThresh1:
    maxFlag = 1
    iSig = smFrame['rank'].argmax() #index of most extreme
    #find indices of extreme
    iSigs = smFrame[smFrame['rank'] > maxThresh1].index
    # connDirect = 'pos'
else:
    maxFlag = 0

if minFlag and maxFlag:
    print 'inconsistant directionality'
    continue

# check if two adjacent concentrations pass thresh2
for iSig in iSigs:
    if maxFlag == 1: 
        threshBool = smFrame['rank'] > maxThresh2
    if minFlag == 1: 
        threshBool = smFrame['rank'] < minThresh2
    iInt = smFrame.index.get_loc(iSig)
    nSm = len(smFrame)
    doseResponse = 0
    #check above and bellow
    if iInt+1 <= nSm-1 and iInt-1 >= 0:
        if sum(threshBool[iInt-1:iInt+2]) == 3:
            doseResponse = 1
    #check two above
    if iInt+2 <= nSm-1:
        if sum(threshBool[iInt:iInt+3]) == 3:
            doseResponse = 1
    #check two bellow
    if iInt-2 >= 0:
        if sum(threshBool[iInt-2:iInt+1]) == 3:
            doseResponse = 1
    if doseResponse == 1:
        continue

if doseResponse == 1:
    print brd + ' dose connection to' + gene1 + ' in ' + cell_tp


#maybe the KD should be the starting query


#sudo code:
#1) if at least one threshold is pass thresh1, continue
#2) if at least two pass thresh2, continue
#3) for the instance that passed thresh1, are there two
# adjacent values that pass thresh2?


