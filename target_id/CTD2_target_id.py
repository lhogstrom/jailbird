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
dg.get_sig_ids(genomic_pert='KD',is_gold=True)
dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')
# dg.test_known_connections(gp_type='KD',metric='spearman',pDescDict=pDescDict,make_graphs=True)
# dg.FDR_correction(pDescDict=pDescDict,metric='spearman',outName='apriori_connections_pass_FDR',alpha=0.2,make_graphs=False)
dg.test_known_connections(gp_type='KD',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='two_sided_dg_graphs_n10M',
                        make_graphs=False,
                        n_rand=10000000,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='KD',
                metric='spearman',
                outName='apriori_two_sided_pass_FDR_n10M',
                alpha=0.2,
                make_graphs=True)
dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR_n10M')
dg.store_parameters_rpt()
# # dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # # dg.test_unknown_rank_product(gp_type='KD')
# # # # dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')
outF = os.path.join(dg.outputdir,'drug-target_summary.txt')
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
                                n_uncorrected=20)
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
                make_graphs=True)
dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR_n10M')
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
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_OE_connection')
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

# # #hyperlink cmd
# outpath = '/'.join([dg.outputdir,'AURKA_gene_to_drug_connections'])
# hyperLnkPath = '/xchip/cogs/web/icmap/hogstrom/AURKA_KD_to_CTD2_connections'
# cmd = ' '.join(['ln -s',
# 		 outpath,
# 		 hyperLnkPath])

params = [] #list of tuples
params.append(('metric',dg.metric))
params.append(('connection_test',dg.connection_test))
params.append(('genomic_perturbation_type',str(dg.gp_type)))
params.append(('is_gold',str(dg.is_gold)))
params.append(('n_permutation',str(dg.n_perm)))
#affogato version
#input target dictionary
# which compounds were actually tested?
# for which comparisons was there not data available?
to.make_rpt(dg.outputdir,params) 



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
            # 'n_connected_cells':self.n_connectedCells}
pFrame = pd.DataFrame(tmpDict)
pFrame['p_log10'] = np.log10(pFrame['p_val'])
pFrame['p_log10_sign'] = pFrame['p_log10'] * pFrame['p_sign']
pFrame = pFrame.sort(columns='p_log10_sign')
plt.plot(pFrame['p_log10_sign'],'.')
plt.title(testGene + ' ' + gp_type)
plt.ylabel('log(p-value)')
plt.xlabel('connection to compounds')
plt.show()
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
    csR = dg.dfCS.ix[brd][ind]
    if gp_type == "OE":
        csR = csR.unstack()
    csR = csR[csR.notnull()]
    rnkR = dg.dfRank.ix[brd][ind]
    if gp_type == "OE":
        rnkR = rnkR.unstack()
    rnkR = rnkR[rnkR.notnull()]
    outF = os.path.join(graphDir,brd +'_' + ind + '_drug-target_summary.txt')
    dg.__make_CS_summary(brd,pDescDict[brd],rnkR,csR,GTDpDict[conn],outF,gp_type,metric)
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
        dg.make_drug_gene_pages(brd,ind,subName,dir_loc=graphName)




