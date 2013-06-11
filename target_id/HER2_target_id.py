#! /usr/bin/env python
'''
analyze drug connections with erbb2 knockdown
'''

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd

work_dir = '/xchip/cogs/projects/target_id/ERBB2_11June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

# Query instances of cps
# CM = mu.CMapMongo()
# pert_List = CM.find({'pert_iname':{'$regex':'5240'}},{'sig_id':True,'pert_iname':True,'pert_id':True,})
# erbb2Lst = CM.find({'pert_iname':{'$regex':'ERBB2'},'pert_type':'trt_sh.cgs'},{'sig_id':True,'pert_iname':True,'pert_id':True,})

#cps to check:
# AZD-8055 - BRD-K69932463
# AS-605240 - BRD-K41895714
# ASG05240
targetDict = {}
targetDict['BRD-K69932463'] = ['ERBB2']
targetDict['BRD-K41895714'] = ['ERBB2']

pDescDict = {}
pDescDict['BRD-K69932463'] = 'AZD-8055'
pDescDict['BRD-K41895714'] = 'AS-605240'

test1 = 'OEB001_A375_96H:BRDN0000399163:-666' #set random sig_id to initialize dgo object
test2 = 'OEB001_A375_96H:BRDN0000400484:-666'

## test OE
reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_OE_connection')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='OE')
# dg.run_drug_gene_query(max_processes=10)
#wait until queries finish
dg.make_result_frames(gp_type='OE')
dg.test_known_connections(pDescDict=pDescDict)
dg.FDR_correction(pDescDict=pDescDict)

### test KD
reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_KD_connection')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD')
# dg.run_drug_gene_query(max_processes=10)
#wait until queries finish
dg.make_result_frames(gp_type='KD')
dg.test_known_connections(gp_type='KD',pDescDict=pDescDict)
dg.FDR_correction(pDescDict=pDescDict)

## # do an umbiased search to see which drugs in CMAP best connect to HER2 CGS


## make summary file
csSer = cpRes[ind]
csSer = csSer[csSer.notnull()]
rnkSer = cpRank[ind]
rnkSer = rnkSer[rnkSer.notnull()]


if not list(rnkSer.index) == list(csSer.index):
    print 'cs and rank series don\'t match'

outF = os.path.join(graphDir + '_drug-target_FDR_summary.txt')    
    self.make(brd,cell,rnkSer,csSer,outF)

outF = os.path.join(graphDir,brd +'_' + ind + '_drug-target_summary.txt')
cp = brd
outfile = outF
gp_type='KD'
    def __make_CS_summary(self,cp,cp_desc,rnkSer,csSer,outfile,gp_type='KD'):
        '''
        make a summary file a given drug-gene pair across cell lines
        inputs:
        rnkSer = pandas series of percent ranks for connections
        csSer = pandas series of cs scores for connections
        '''
        if gp_type == 'KD':
            headers = ['query_sig','cp_pert_iname','CGS_KD','cell','cs', 'RnkPerc']
        if gp_type == 'OE':
            headers = ['query_sig','cp_pert_iname','OE_compared','cell','cs', 'percent_rank']
        with open(outfile,'w') as f:
            f.write('\t'.join(headers) + '\n')
            for i,cs in enumerate(csSer):
                query = rnkSer.index[i]
                desc = pDescDict[cp]
                rank = rnkSer[i]
                #query, target gene, cs, rank, RnkPerc
                f.write('\t'.join([query,cp_desc,rnkSer.name,query.split('_')[1],str(cs),str(rank)[:7]]) + '\n')

        outF = os.path.join(celldir,cell1 + '_drug-target_FDR_summary.txt')
        
        with open(outF,'w') as f:
            f.write('\t'.join(headers) + '\n')
            for pert in targetRnkPercs[cell1]:
                qInds = queryInd[pert]
                for target in targetRnkPercs[cell1][pert]:
                    RnkPercs = targetRnkPercs[cell1][pert][target]
                    for i,RnkPerc in enumerate(RnkPercs):
                        if RnkPerc <= pMaxThresh:
                            query = resCids[qInds[i]]
                            cs = targetCS[cell1][pert][target][i]
                            rank = targetRanks[cell1][pert][target][i]
                            #query, target gene, cs, rank, RnkPerc
                            f.write('\t'.join([query,pDescDict[pert],target,str(cs),str(rank),str(RnkPerc)[:7]]) + '\n')


        #set output dir
        gp_type='KD'
        n_rand=100000
        make_graphs=True
        graphDir = dg.outputdir + '/drug_target_graphs'
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
        cols = dg.dfRank.columns
        pDict = {}
        pVec = []
        prog = progress.DeterminateProgressBar('Connection test')
        for ibrd,brd in enumerate(dg.targetDict):
            # skip pert if not in result file
            prog.update(brd,ibrd,len(dg.targetDict))
            if not brd in brdRsltSet:
                brdSkipped.append(brd)
                continue
            targets = dg.targetDict[brd]
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
                    rnkSer = cpRank[ind]/100
                    rnkSer = rnkSer[rnkSer.notnull()]
                    #skip if cgs not tested in the same cell line as cp
                    if len(rnkSer) == 0:
                        cgsSkipped.append(ind)
                        continue
                    ### calculate p-value - based on percent rank products
                    testStat = rnkSer.prod()
                    n_obs = rnkSer.shape[0]
                    # theoretical null
                    ### simulate random draws from percent rank list
                    permMtrx = np.random.rand(n_obs,n_rand)
                    nullDist = permMtrx.prod(axis=0)
                    #number of null values more extreme than observed (one sided)
                    exVals = nullDist[nullDist<testStat]
                    nExtreme = len(exVals)
                    pVal = (nExtreme+1)/float(len(nullDist))
                    pVec.append(pVal)
                    pDict[brd + '-' + ind] = pVal
                    # if valCounts[ind] > 9:
                    # if pVal < .01:
                    if make_graphs:
                        ### cs wadden gram
                        sKeysStr = []
                        count = 0
                        for i,cs in enumerate(cpRes[ind]):
                            if pd.isnull(cs):
                                continue
                            else:
                                count = count + 1
                                sKeysStr.append(cpRes.index[i].split('_')[1])
                                yVals = count
                                plt.scatter(cs,yVals)
                        plt.xlim((-1, 1))
                        plt.ylim((0,count+1))
                        plt.yticks(range(1, count + 2), sKeysStr, rotation = 0)
                        plt.xlabel('wtcs')
                        plt.ylabel('cell line')
                        plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type)
                        plt.savefig(os.path.join(graphDir,brd +'_' + ind + '_connections.png'))
                        plt.close()
                        #rank wadden gram
                        sKeysStr = []
                        count = 0
                        rnkList = cpRank[ind]
                        for i,rnk in enumerate(cpRank[ind]):
                            if pd.isnull(rnk):
                                continue
                            else:
                                count = count + 1
                                sKeysStr.append(cpRes.index[i].split('_')[1])
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
        self.brdSkipped = brdSkipped
        self.cgsSkipped = cgsSkipped
        self.pDict = pDict
        self.pVec = pVec



### trouble shoot pd.concat issue
        gp_type='OE'
        work_dir2 = work_dir + '/drug_OE_connection'
        cellDirs = [f for f in os.listdir(work_dir2) if os.path.isdir(work_dir2+'/'+f)]
        prog = progress.DeterminateProgressBar('dataframe read')
        df = pd.DataFrame()
        dfRank = pd.DataFrame()
        #loop through each cell line add to df
        for icell, cell1 in enumerate(cellDirs):
            #define directories and load in outputs
            outdir = os.path.join(work_dir2,cell1,'sig_query_out')
            if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx'):
                print cell1 + 'no query result file'
                continue #if no results file, skip loop
            rsltFile = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
            rslt = gct.GCT()
            rslt.read(rsltFile)
            prog.update(cell1,icell,len(cellDirs))
            rsltF = rslt.frame
            rsltF = rsltF.T
            indVals = rsltF.index.values
            pertVals = [ind.split(':')[1][:13] for ind in indVals]
            #make the column name gene and pert time
            geneVals = []
            for ind in rsltF.columns:
                if gp_type == 'KD':
                    gene = ind.split(':')[1]
                if gp_type == 'OE':
                    brdn = ind.split(':')[1]
                    gene = dg.BRDNdict[brdn]
                tp = ind.split(':')[0].split('_')[-1]
                gname = '_'.join([gene, tp])
                geneVals.append(gname)
            if len(geneVals) > len(set(geneVals)):
                print 'duplicate CGS for this celline'
            newF = rsltF
            newF.index = [pertVals, rsltF.index.values]
            newF.columns = geneVals
            rankF = newF.rank(ascending=False,axis=1)
            perRankF = rankF / float(rankF.shape[1]) * 100.0
            #add cell line result to combined df
            if len(df) == 0:
                df = newF
                dfRank = perRankF
            else:
                df = pd.concat([df,newF],axis=0,join='outer')
                dfRank = pd.concat([dfRank,perRankF],axis=0)

df1 = df[:1]
newF1 = newF[:1]

df1 = df.ix['BRD-K41895714']
newF1 = newF.ix['BRD-K41895714']
df2 = pd.concat([df1,newF1],axis=0)

item = 'A2M_96H'
print [x for x in newF.columns if x == item]
import collections
print [x for x, y in collections.Counter(list(newF.columns)).items() if y > 1]


