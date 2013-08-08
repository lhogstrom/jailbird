#! /usr/bin/env python
'''
summly results for groups of compounds from ctd2_annots

make heatmaps of summly rank for each class of compounds 

Author : Larson Hogstrom
Date: Aug 2013

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
import matplotlib

work_dir = '/xchip/cogs/hogstrom/analysis/summly/cp_class'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

drugFile = '/xchip/cogs/projects/target_id/ctd2_annots/ctd2_merged_mapped_genes.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
targetSet = set(drugLabels['gene_dw_annot'])
### pair brd and pert_iname
allCtd2 = drugLabels['pert_id']
allinames = drugLabels['pert_iname']
pDescDict = {}
for ibrd,brd in enumerate(allCtd2):
    pDescDict[brd] = allinames[ibrd]

### run queries for cps
# cgsCells = ['A375', 'A549', 'ASC', 'HA1E', 'HEPG2', 'HCC515', 'HT29', 'MCF7', 'NPC', 'PC3', 'VCAP']
# processes = set()
# file1 = work_dir + '/' + 'ctd2_sigs.grp'
# for pert in allCtd2:
#     sigs = []
#     for cell in cgsCells:
#         CM = mu.CMapMongo()
#         #dose min
#         # pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell,'pert_dose':{'$gt':5}},{'sig_id':True},limit=1)
#         #no dose min
#         pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell},{'sig_id':True},limit=1)
#         if pert_query:
#             sigs.append(pert_query[0])
#         # else:
#         #     print 'no sig for ' + pert + ' in ' + cell
#     with open(file1, 'a') as f:
#         [f.write(x + '\n') for x in sigs]
# metric = 'wtcs'
# queryDir = work_dir + '/ctd2_sig_query'
# if not os.path.exists(queryDir):
#     os.mkdir(queryDir)
# cmd = ' '.join(['rum -q local -f sig_query_tool',
#          '--sig_id ' + file1,
#          '--metric ' + metric,
#          '--column_space full',
#          '--out ' + queryDir,
#          '--mkdir false',
#          '--save_tail false'])
#          # '--row_space bing', 
# os.system(cmd)

########################################
###### full ctd2 summly results ######
########################################

brdGrpList = []
grpSet = set(drugLabels['gene_dw_annot'])
grpToCp = {}
for grp in grpSet:
    grpPerts = drugLabels['pert_id'][drugLabels['gene_dw_annot'] == grp]
    grpToCp[grp] = list(grpPerts.values)
    brdGrpList.extend(grpPerts.values)
# get list of cps in summly dir
basePath = work_dir + '/ctd2_sig_query'
# basePath = '/xchip/cogs/hogstrom/analyszis/summly/cp_class/ctd2_sig_query'
dateID = 'aug01/my_analysis.sig_summly_tool.2013080119394091'
summDir = '/'.join([basePath,dateID])
cpDirs = [f for f in os.listdir(summDir) if os.path.isdir(summDir+'/'+f)]

graphDir = work_dir + '/ctd2_graph_out'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
### examine each functional gorup
for grpGene in grpToCp:
    if grpGene == '-666':
        continue
    grp = grpToCp[grpGene]
    if not grp: # skip if grp is empty
        continue
    # matrices for grp connections
    nGrp = len(grp)
    grp_sum_score = np.zeros((nGrp,nGrp))
    grp_PercSummly = np.zeros((nGrp,nGrp))
    grp_rank = np.zeros((nGrp,nGrp))
    for ibrd,brd in enumerate(grp):
        # brd = 'BRD-A15100685'
        basePath = work_dir + '/sig_query_10um'
        dateID = 'jul29/my_analysis.sig_summly_tool.2013072914520591'
        inFile = '/'.join([summDir,
                        brd,
                        brd+'_summly.txt'])
        sumRes = pd.io.parsers.read_csv(inFile,sep='\t')
        # filter to only cps / cgs
        pd.io.parsers.read_csv
        cpRes = sumRes[sumRes['pert_type'] == 'trt_cp']
        cpRes['rank'] = np.arange(1,len(cpRes)+1)
        cgsRes = sumRes[sumRes['pert_type'] == 'trt_sh.cgs']
        cgsRes['rank'] = np.arange(1,len(cgsRes)+1)
        oeRes = sumRes[sumRes['pert_type'] == 'trt_oe']
        oeRes['rank'] = np.arange(1,len(oeRes)+1)
        #check group connection 
        for ibrd2, brd2 in enumerate(grp):
            indSum = cpRes[cpRes['pert_id'] == brd2]['sum_score']
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
    ### print group heatmap
    fig = plt.figure(1, figsize=(20, 8))
    plt.suptitle(grpGene + ' compound group',fontsize=14, fontweight='bold')
    plt.subplot(121)
    plt.title('percent summly rank')
    plt.imshow(grp_PercSummly,
            interpolation='nearest',
            cmap=matplotlib.cm.RdBu_r,
            vmin=0, 
            vmax=1)
    ytcks = [pDescDict[x] for x in grp]
    plt.xticks(np.arange(len(grp)), ytcks,rotation=75)
    plt.yticks(np.arange(len(grp)),ytcks)
    plt.colorbar()
    plt.subplot(122)
    plt.title('sum_score')
    plt.imshow(grp_sum_score,
            interpolation='nearest',
            cmap=matplotlib.cm.RdBu_r,
            vmin=-1, 
            vmax=1)
    plt.xticks(np.arange(len(grp)), ytcks,rotation=75)
    # ytcks = [pDescDict[x] for x in avicinsBrds]
    plt.yticks(np.arange(len(grp)),ytcks)
    plt.colorbar()
    outF = os.path.join(graphDir,grpGene + '_compound_group_heatmap.png')
    fig.savefig(outF, bbox_inches='tight')
    plt.close()



### average number of hairpins that go into a CGS at the top of the summmly lists
genes = list(cgsRes['pert_iname'].values)
CM = mu.CMapMongo()
nrepFrame = CM.find({'pert_iname' : {'$in' : genes},'pert_type':'trt_sh.cgs'}, {'pert_iname':True,'cell_id':True,'distil_nrep':True}, toDataFrame = True)
nrepFrame.sort_index(by = 'pert_iname', inplace = True)
#calculate average number of distil_nreps for each cgs across cell lines
avFrame = nrepFrame.groupby(['pert_iname']).mean()
avNrep = avFrame['distil_nrep']

upReps = []
dnReps = []
for ibrd,brd in enumerate(cpDirs):
    basePath = work_dir + '/sig_query_10um'
    dateID = 'jul29/my_analysis.sig_summly_tool.2013072914520591'
    inFile = '/'.join([summDir,
                    brd,
                    brd+'_summly.txt'])
    sumRes = pd.io.parsers.read_csv(inFile,sep='\t')
    # filter to only cps / cgs
    pd.io.parsers.read_csv
    cpRes = sumRes[sumRes['pert_type'] == 'trt_cp']
    cpRes['rank'] = np.arange(1,len(cpRes)+1)
    cgsRes = sumRes[sumRes['pert_type'] == 'trt_sh.cgs']
    cgsRes['rank'] = np.arange(1,len(cgsRes)+1)
    oeRes = sumRes[sumRes['pert_type'] == 'trt_oe']
    oeRes['rank'] = np.arange(1,len(oeRes)+1)
    #sort by summly results
    resNrep = avNrep.reindex(index=cgsRes['pert_iname'])
    up50 = resNrep.ix[:25].values
    dn50 = resNrep.ix[-25:].values
    upReps = np.append(up50,upReps)
    dnReps = np.append(dn50,dnReps)

### make histogram of nrep counts - baseline vs summly results
fig = plt.figure(1, figsize=(20, 8))
plt.suptitle('distil_nrep distribution',fontsize=14, fontweight='bold')
plt.subplot(211)
plt.title('all KDs')
plt.xlabel('average distil_nreps per gene')
n, bins, patches = plt.hist(avNrep.values,40)
plt.subplot(212)
plt.title('CTD2 summly top connecting KDs')
nR, binsR, patchesR = plt.hist(upReps,bins=bins,color='r')
plt.xlabel('average distil_nreps per gene')

#plot ratio
plt.plot(bins[:40],nrepRatio,'o')
plt.title('overrepresentation of high distil_nreps in summly results',fontsize=14, fontweight='bold')
plt.xlabel('average number of hairpins in cgs (distil_nreps)')
plt.ylabel('ratio top connectors to baseline')





### scratch
# inameSer = nrepFrame['pert_iname']
# nrepFrame['distil_nrep']
# inameSer.value_counts()
# akt1 = nrepFrame[nrepFrame['pert_iname'] == 'AKT1']

# scratchFrm = pd.Series(list(set(nrepFrame['pert_iname'].values)))
# scratchFrm = scratchFrm.sort()
# outF = '/xchip/cogs/hogstrom/analysis/summly/cp_class/KD_gene_symbols.txt'
# scratchFrm.to_csv(path=outF,sep='\t')
# outF = gseaDir + '/' + brd + '_KD_summly_result.txt'
# scoreSeries = pd.Series(cgsRes['sum_score'].values, index=cgsRes['pert_iname'])
# scoreSeries.to_csv(path=outF,sep='\t')

# n, bins, patches = plt.hist(avNrep.values,40)
# plt.hist(up50,bins=bins)

