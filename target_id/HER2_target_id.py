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

work_dir = '/xchip/cogs/projects/target_id/ERBB2_12June2013'
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
# reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_OE_connection')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='OE')
dg.run_drug_gene_query(max_processes=10)
#wait until queries finish
dg.make_result_frames(gp_type='OE')
dg.test_known_connections(pDescDict=pDescDict)
dg.FDR_correction(pDescDict=pDescDict)

### test KD
# reload(dgo)
# dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_KD_connection')
# dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD')
# dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
# dg.make_result_frames(gp_type='KD')
# dg.test_known_connections(gp_type='KD',pDescDict=pDescDict)
# dg.FDR_correction(pDescDict=pDescDict)

## # do an umbiased search to see which drugs in CMAP best connect to HER2 CGS




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
                df = pd.concat([df,newF],axis=0)
                dfRank = pd.concat([dfRank,perRankF],axis=0)

### 
h1 = df
h1.index.names = ['cp','sig_id']
p1 = newF
p1.index.names = ['cp','sig_id']
pc3Res = '/xchip/cogs/projects/target_id/ERBB2_12June2013/drug_OE_connection/PC3/sig_query_out/result_WTCS.LM.COMBINED_n4x2937.gctx'
ha1eRes = '/xchip/cogs/projects/target_id/ERBB2_12June2013/drug_OE_connection/HA1E/sig_query_out/result_WTCS.LM.COMBINED_n4x2159.gctx'
prslt = gct.GCT()
prslt.read(pc3Res)
hrslt = gct.GCT()
hrslt.read(ha1eRes)
zh = h1.ZRSR2_96H
zp = p1.ZRSR2_96H
ah = h1.ABAT_96H
ap = p1.ABAT_96H

g = pd.concat([zh,ph])
f = pd.concat([ap,zp])

zc = pd.concat([zh,ph],axis=0)

#re-index series
cpBRDs = zh.index.levels[0]

# hF = pd.concat([ah,zh],axis=1)
# pF = pd.concat([ap,zp],axis=1)
hF = h1.ix[:,:3]
pF = p1.ix[:,:3]
Cnct = pd.concat([hF,pF],axis=0)


df1 = df[1:3]
newF1 = newF[1:3]

df1 = df.ix['BRD-K41895714']
newF1 = newF.ix['BRD-K41895714']
df2 = pd.concat([df1,newF1],axis=0)

item = 'A2M_96H'
print [x for x in newF.columns if x == item]
import collections
print [x for x, y in collections.Counter(list(newF.columns)).items() if y > 1]

#toy df
dF = pd.DataFrame(np.random.randn(10, 4))
pieces = [dF[:3], dF[3:7], dF[7:]]
concatenated = pd.concat(pieces, keys=['first', 'second', 'third'])

pieces = [dF.ix[:, [0, 1]], dF.ix[:, [2]], dF.ix[:, [3]]]
result = pd.concat(pieces, axis=1, keys=['o','t','th'])

