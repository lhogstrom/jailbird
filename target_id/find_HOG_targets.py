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

work_dir = '/xchip/cogs/projects/target_id/HOG_8July2013_not_gold'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

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
### make query with only HOG plates
####################################
dg = dgo.QueryTargetAnalysis(out=work_dir)
pert_list = list(hogSet)
is_gold =False
genomic_pert = 'KD'
brdCounts = []
fullPertList = []
### for each drug perturbations of interest - find all instances in CMAP
prog = progress.DeterminateProgressBar('perturbation cid query')
if pert_list:
    for i,pert in enumerate(pert_list):
        prog.update('querying cps',i,len(pert_list))
        CM = mu.CMapMongo()
        if is_gold == True:
            pert_query = CM.find({'sig_id':{'$regex':'HOG'},'pert_id':{'$regex':pert},'is_gold':True},{'sig_id':True,'cell_id':True})
        if is_gold == False:
            pert_query = CM.find({'sig_id':{'$regex':'HOG'},'pert_id':{'$regex':pert}},{'sig_id':True,'cell_id':True})
        if pert_query:
            brdCounts.append(len(pert_query))
            fullPertList.extend(pert_query)
cell_lines_tested = []
cellsAll = [sig['cell_id'] for sig in fullPertList]
uniqCells = list(set(cellsAll))
prog = progress.DeterminateProgressBar('genomic pert query')
### 1) which celll lines have tested the target with a genomic pert - write the cids to a file
### 2) write cp cig ids to a file if there are CGSs in that cell line
if genomic_pert == 'OE':
    BRDNdict = {} #dictionary of gene symbol for each BRDN number
for i,cell1 in enumerate(uniqCells):
    #get all CGS for a cell line
    prog.update('querying',i,len(uniqCells))
    CM = mu.CMapMongo()
    if genomic_pert == 'KD':
        CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','is_gold':True,'cell_id':cell1},{'sig_id':True,'pert_iname':True})
    if genomic_pert == 'OE':
        CGSbyCell = CM.find({'pert_type':'trt_oe','is_gold':True,'cell_id':cell1},{'sig_id':True,'pert_iname':True})
        if CGSbyCell:
            for query in CGSbyCell:
                brdn = query['sig_id'].split(':')[1]
                BRDNdict[brdn] = query['pert_iname']
    if CGSbyCell:
        cell_lines_tested.append(cell1)
        outdir = os.path.join(work_dir,cell1)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        nCGS = len(CGSbyCell)
        sigF = os.path.join(outdir, cell1+ '_genomic_sig_ids_n' + str(nCGS) + '.grp')
        with open(sigF, 'w') as f:
            for sig in CGSbyCell:
                f.write(sig['sig_id'] + '\n')
        sigIDlist = []
        for sig in fullPertList:
            if sig['cell_id'] == cell1:
                sigIDlist.append(sig['sig_id'])
        sigIDlist = list(set(sigIDlist))
        #write drug signatures by cell line to a file
        sigF = os.path.join(outdir,cell1 + '_cp_sig_ids.grp')
        with open(sigF, 'w') as f:
            [f.write(x + '\n') for x in sigIDlist]
dg.cell_lines_tested = cell_lines_tested
if genomic_pert == 'OE':
    dg.BRDNdict = BRDNdict

  
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
                        n_rand=1000000,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='KD',
                metric='spearman',
                outName='apriori_FDR_pass',
                alpha=0.2,
                make_graphs=True,
                specificity_graph=True)
dg.fdr_html_summary(fdrDir='test_FDR2',specificity_graph=True)


# reload module
# dgCopy = dg
reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir)
dg.add_dictionary(targetDict=hogTargetDict)
dg.dfRank = dgCopy.dfRank
dg.dfCS = dgCopy.dfCS
####################################
### make dose specific df ##########
####################################
sigIDs = [x[1] for x in dg.dfRank.index]
#make a pandas series of sigIDs and doses
doses = [float(x.split(':')[-1]) for x in sigIDs]
doseSer = pd.Series(doses,index=sigIDs)
#find the most frequent doses
doseCount = doseSer.value_counts()[:9]
doseGroup = doseCount.index
doseGroup = doseGroup.values
doseGroup.sort()
doseGroup = list(doseGroup)
#round the each dose to one of the 9 most frequent values
nearest_doses = []
for dose in doseSer:
	doseNear = doseGroup[(np.abs(doseGroup-dose)).argmin()]
	nearest_doses.append(doseNear)
#frame of origonal and rounded doses
doseFrame = pd.DataFrame({'dose': doseSer,'dose_round':nearest_doses})
# #build a smaller frame of doses
# for dose in doseGroup:
# 	smFrame = doseFrame[doseFrame['dose_round'] == dose]

doseFrameRnk = dg.dfRank
doseFrameRnk['dose'] = doseSer.values
doseFrameRnk['dose_round'] = nearest_doses
doseFrameCS = dg.dfCS
doseFrameCS['dose'] = doseSer.values
doseFrameCS['dose_round'] = nearest_doses

# calculate the number of connections over spearman .4
metricthresh = .5
doseCounts = []
threshCounts = []
relativeThrCnts = []
for dose in doseGroup:
	smFrameCS = doseFrameCS[doseFrameCS['dose_round'] == dose]
	doseCount = smFrameCS.shape[0]
	doseCounts.append(doseCount)
	threshFrame = smFrameCS[smFrameCS >= metricthresh] 
	threshCount = sum(threshFrame.count())
	threshCounts.append(threshCount)
	relativeThrCnt = threshCount/float(doseCount*smFrameCS.shape[1])
	relativeThrCnts.append(relativeThrCnt)

# graph 
plt.bar(np.arange(len(relativeThrCnts)),relativeThrCnts)
plt.title('percent of all drug-CGS pairwise comparisons\n passing .5 spearman threshold')
plt.ylabel('percent of all pairwise comparisons')
plt.xlabel('dose - um')
doseStrs = [str(x) for x in doseGroup]
plt.xticks(np.arange(len(relativeThrCnts))+.5, doseStrs, rotation = 45)

# is_gold
for kd_tp in dg.dfRank.columns:
	grouped = doseFrameCS[kd_tp].groupby(doseFrameCS['dose_round'])
	grpMeans = grouped.mean()
	plt.plot(grpMeans.values,alpha=.1)
plt.ylabel('percent of all pairwise comparisons')
plt.xlabel('dose - um')
doseStrs = [str(x) for x in doseGroup]
plt.xticks(np.arange(len(relativeThrCnts))+.5, doseStrs, rotation = 45)
plt.show()

### check dose paterns of expected connections 

#jame st