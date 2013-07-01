 #! /usr/bin/env python

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd

work_dir = '/xchip/cogs/projects/target_id/CTD2_25June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

### make target_dict
# targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_short.txt'
targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_drug_targets.txt'
targetDict = {}
pDescDict = {}
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

### test KD
reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')
# dg.test_known_connections(gp_type='KD',metric='spearman',pDescDict=pDescDict,make_graphs=True)
# dg.FDR_correction(pDescDict=pDescDict,metric='spearman',outName='apriori_connections_pass_FDR',alpha=0.2,make_graphs=False)
dg.test_known_connections(gp_type='KD',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='test_dg_graphs',
                        conn_thresh=.05,
                        make_graphs=True,
                        n_rand=100000,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='KD',
                metric='spearman',
                outName='test_FDR2',
                alpha=0.2,
                make_graphs=True,
                specificity_graph=True)
dg.fdr_html_summary(fdrDir='test_FDR2',specificity_graph=True)
# dg.store_parameters_rpt()
# # dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # # dg.test_unknown_rank_product(gp_type='KD')
# # # # dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')
outF = os.path.join(dg.outputdir,'drug-target_summary_peyton.txt')
dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')

# #random
# geneList = ['ERBB2','MUC1','PIK3CA','MTOR','PPARG']
# for gene in geneList:
# 	dg.gene_to_drug_similarity(testGene=gene,gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# #AURKA
# geneList = ['AURKA','AURKB','AURKAIP1']
inFile = '/xchip/cogs/projects/target_id/CTD2_25June2013/genes_with_connections.txt'
cgsList = []
with open(inFile,'rt') as f:
    for string in f:
        splt = string[:-1]
        cgsList.append(splt)
geneList = set(cgsList)
for gene in geneList:
    dg.gene_to_drug_similarity(testGene=gene,
                                gp_type='KD',
                                metric='spearman',
                                outName='gene_to_drug_connections',
                                pDescDict=pDescDict,
                                n_rand=10000,
                                n_uncorrected=20,
                                connection_test='two_sided')
# #drugbank connections
# geneList = ['PPARG,', 'FKBP1A', 'KIF11', 'MTOR', 'HMGCR', 'RRM1', 'ESR1', 'NR3C1', 'HMGCR', 'NNR3C1', 'HMGCR', 'NR3C1', 'PSMB1', 'PSMB5', 'RAF1', 'BRAF', 'CDK4', 'ESR1', 'NR3C1', 'NR3C1', 'NR3C1', 'NR3C1', 'R3C1', 'EGFR', 'HMGCR','EGFR', 'RRM1']


# # ### TEST OE
# # reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_OE_connection')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='OE',is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='OE',metric='spearman')
dg.test_known_connections(gp_type='OE',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='two_sided_dg_graphs',
                        n_rand=10000000,
                        make_graphs=False,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='OE',
                metric='spearman',
                outName='apriori_two_sided_pass_FDR_n10M',
                alpha=0.2,
                make_graphs=True,
                specificity_graph=True)
dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR_n10M',specificity_graph=True)
# dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # dg.test_unknown_rank_product(gp_type='KD')
# # # # dg.FDR_correction(pDe
outF = os.path.join(dg.outputdir,'drug-target_two_sided_summary.txt')
dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')
dg.store_parameters_rpt()


# # #re-asign varabile so computation is not lost when reloading the class
# # #dgCopy = dg
reload(dgo)
# dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
dg.dfRank = dgCopy.dfRank
dg.dfCS = dgCopy.dfCS
dg.pVec = dgCopy.pVec
dg.pDict = dgCopy.pDict
dg.countDict = dgCopy.countDict
dg.cellCountDict = dgCopy.cellCountDict
dg.nConnectionDict = dgCopy.nConnectionDict
dg.connectionsPassFDR = dgCopy.connectionsPassFDR
dg.pThreshFDR = dgCopy.pThreshFDR
dg.connection_test = dgCopy.connection_test

# # #hyperlink cmd
# outpath = '/'.join([dg.outputdir,'AURKA_gene_to_drug_connections'])
# hyperLnkPath = '/xchip/cogs/web/icmap/hogstrom/AURKA_KD_to_CTD2_connections'
# cmd = ' '.join(['ln -s',
# 		 outpath,
# 		 hyperLnkPath])


### make specificity graph
testGene='MTOR'
gp_type='KD'
metric='spearman'
pDescDict=pDescDict
outName='gene_to_drug_connections'
n_rand=100000
make_graphs=True
connection_test = 'two_sided'

graphName = testGene + '_' + outName
graphDir = dg.outputdir + '/' + testGene + '_' + outName
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
#find genomic perts that match input
geneSymMtch = []
igeneSymMtch = []
for icol, col in enumerate(dg.dfRank.columns):
    if col.split('_')[0] == testGene:
        geneSymMtch.append(col)
        igeneSymMtch.append(icol)
    else:
        continue
GTDpVec = []
GTDpDict = {}
GTDsignDict = {}
# identify percent ranks of each drug-gene combo
for ind in geneSymMtch:
    geneRnkSer = dg.dfRank[ind]
    geneRnkSer = geneRnkSer[geneRnkSer.notnull()]
    geneCSSer = dg.dfCS[ind]
    geneCSSer = geneCSSer[geneCSSer.notnull()]
    cpList = [g[0] for g in geneRnkSer.index]
    cpSet = list(set(cpList))
    for brd in cpSet:
        rnkSer = geneRnkSer[brd]
        csSer = geneCSSer[brd]
        ### calculate p-value - based on percent rank products
        rnkSm = rnkSer/100
        testStat = rnkSm.prod()
        n_obs = rnkSer.shape[0]
        # theoretical null
        ### simulate random draws from percent rank list
        permMtrx = np.random.rand(n_obs,n_rand)
        nullDist = permMtrx.prod(axis=0)
        #number of null values more extreme than observed (one sided)
        exVals = nullDist[nullDist<testStat]
        nExtreme = len(exVals)
        pVal = (nExtreme+1)/float(len(nullDist))
        ### one/ two sided p-value
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
                pSign = 1
            else:
                pVal = (NnExtreme+1)/(float(len(nullDist))/2)
                pSign = -1
        #save p-values to dictionary
        GTDpVec.append(pVal)
        GTDpDict[brd + ':' + ind] = pVal
        GTDsignDict[brd + ':' + ind] = pSign
dg.GTD_pDict = GTDpDict
dg.GTD_testGene = testGene 
#create a pandas series and order p-values
pSer = pd.Series(GTDpDict)
pSer.sort()
tmpDict = {'p_val':GTDpDict,
            'p_sign':GTDsignDict}
            # 'n_connected_cells':dg.n_connectedCells}
pFrame = pd.DataFrame(tmpDict)
pFrame['p_log10'] = np.log10(pFrame['p_val'])
pFrame['p_log10_sign'] = pFrame['p_log10'] * pFrame['p_sign']
pFrame = pFrame.sort(columns='p_log10_sign')
#make plot
plt.plot(pFrame['p_log10_sign'],'.')
plt.title(testGene + ' ' + gp_type)
plt.ylabel('log(p-value)')
plt.xlabel('connection to compounds')
plt.savefig(os.path.join(graphDir,testGene + '_specificity_rank_product_pValues.png'))
plt.close()
#chose top/ bottom connections to graph
examine = pSer[:n_uncorrected]
examine = examine.append(pSer[-n_uncorrected:])
for conn in examine.index:
    brd = conn.split(':')[0]
    if brd[:3] != 'BRD':
        print 'skipped non-brd field: ' + brd
        continue
    ind = conn.split(':')[-1]
    #copied from FDR mdule:
    ### cs wadden gram
    csR = self.dfCS.ix[brd][ind]
    if gp_type == "OE":
        csR = csR.unstack()
    csR = csR[csR.notnull()]
    rnkR = self.dfRank.ix[brd][ind]
    if gp_type == "OE":
        rnkR = rnkR.unstack()
    rnkR = rnkR[rnkR.notnull()]
    outF = os.path.join(graphDir,brd +'_' + ind + '_drug-target_summary.txt')
    self.__make_CS_summary(brd,pDescDict[brd],rnkR,csR,GTDpDict[conn],outF,gp_type,metric)
    sKeysStr = []
    count = 0
    for i,cs in enumerate(csR):
        if pd.isnull(cs):
            continue
        else:
            count = count + 1
            if gp_type == 'KD':
                sKeysStr.append(csR.index[i].split('_')[1])
            if gp_type == 'OE':
                sKeysStr.append(csR.index[i][0].split('_')[1])
            yVals = count
            plt.scatter(cs,yVals)
    plt.xlim((-1, 1))
    plt.ylim((0,count+1))
    plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
    plt.xlabel(metric)
    plt.ylabel('cell line')
    plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type)
    plt.savefig(os.path.join(graphDir,brd +'_' + ind + '_connections.png'))
    plt.close()
    #rank wadden gram
    sKeysStr = []
    count = 0
    # rnkList = cpRank[ind]
    for i,rnk in enumerate(rnkR):
        if pd.isnull(rnk):
            continue
        else:
            count = count + 1
            if gp_type == 'KD':
                sKeysStr.append(csR.index[i].split('_')[1])
            if gp_type == 'OE':
                sKeysStr.append(csR.index[i][0].split('_')[1])
            yVals = count
            plt.scatter(rnk,yVals)
    plt.xlim((0, 100))
    plt.ylim((0,count+1))
    plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
    plt.xlabel('percent rank')
    plt.ylabel('cell line')
    plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type)
    plt.savefig(os.path.join(graphDir,brd +'_' + ind + '_percent_rank.png'))
    plt.close()
### make html summary page
indexfile = os.path.join(graphDir,'index.html')
with open(indexfile,'w') as f:
    lineWrite = '<h2><CENTER>Conenctions to gene ordered by p-value <CENTER></h2>'
    f.write(lineWrite + '\n')
with open(indexfile,'a') as f:
    for conn in examine.index:
        brd = conn.split(':')[0]
        if brd[:3] != 'BRD':
            print 'skipped non-brd field: ' + brd
            continue
        ind = conn.split(':')[-1]
        prcRnkName = '_'.join([brd,ind,'percent_rank.png'])
        # csF = os.path.join(graphDir,prcRnkName)
        # lineWrite = '<BR>' + '<tr><td><img src=' + csF + '></td></tr>'
        subName = '_'.join([brd,ind +'.html'])
        # subPage = os.path.join(graphDir,subName)
        lineWrite = '<p><a href="' + subName + '"><img src="' + prcRnkName + '"alt="' + brd + '"></a></p>'
        f.write(lineWrite + '\n')
        #make sub-page
        self.make_drug_gene_pages(brd,ind,subName,dir_loc=graphName)



pSer = pd.Series(GTDpDict)
pSer.sort()
tmpDict = {'p_val':GTDpDict,
            'p_sign':GTDsignDict}
            # 'n_connected_cells':self.n_connectedCells}
pFrame = pd.DataFrame(tmpDict)
pFrame['p_log10'] = np.log10(pFrame['p_val'])
pFrame['p_log10_sign'] = pFrame['p_log10'] * pFrame['p_sign']
pFrame = pFrame.sort(columns='p_log10_sign')
#make plot
plt.plot(pFrame['p_log10_sign'],'.')
plt.title(testGene + ' ' + gp_type)
plt.ylabel('log(p-value)')
plt.xlabel('connection to compounds')
# plt.show()
combo = 'BRD-K03618428:MTOR_96H'
pLogSer = pFrame['p_log10_sign']
vcombo = pLogSer[combo]
icombo = pLogSer.index.get_loc(combo)
plt.plot(icombo,vcombo,'r.')
plt.show()
# plt.savefig(os.path.join(graphDir,testGene + '_specificity_rank_product_pValues.png'))
# plt.close()


gp_type='KD'
metric='spearman'
pDescDict=pDescDict
outName='test_dg_graphs'
make_graphs=False
n_rand=100000
connection_test='two_sided'
conn_thresh=.01
### cell dict
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
                connCellDict[brd + ':' + ind] = cells
            nConnectionDict[brd + ':' + ind] = n_connections
            #make summary output
            outF = os.path.join(graphDir,brd +'_' + ind + '_drug-target_summary.txt')
            # dg.__make_CS_summary(brd,pDescDict[brd],rnkSer,csSer,pVal,outF,gp_type,metric)
            if make_graphs:
                ### cs wadden gram
                csGraph = os.path.join(graphDir,brd +'_' + ind + '_connections.png')
                dg.__connection_dot_plot(csSer,
                                        pDescDict[brd],
                                        ind,
                                        csGraph,
                                        gp_type='KD',
                                        axis=metric)
                ### percent rank wadden gram
                csGraph = os.path.join(graphDir,brd +'_' + ind + '_percent_rank.png')
                dg.__connection_dot_plot(rnkSer,
                                        pDescDict[brd],
                                        ind,
                                        csGraph,
                                        gp_type='KD',
                                        axis='perc_rank')

