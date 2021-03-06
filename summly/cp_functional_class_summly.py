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

drugFile = '/xchip/cogs/projects/cp_annot/ctd2_annots/ctd2_merged_mapped_genes.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
targetSet = set(drugLabels['gene_dw_annot'])

allCtd2 = drugLabels['pert_id']
allinames = drugLabels['pert_iname']
pDescDict = {}
for ibrd,brd in enumerate(allCtd2):
    pDescDict[brd] = allinames[ibrd]

### run queries for cps
cgsCells = ['A375', 'A549', 'ASC', 'HA1E', 'HEPG2', 'HCC515', 'HT29', 'MCF7', 'NPC', 'PC3', 'VCAP']
processes = set()
file1 = work_dir + '/' + 'ctd2_sigs.grp'
for pert in allCtd2:
    sigs = []
    for cell in cgsCells:
        CM = mu.CMapMongo()
        #dose min
        # pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell,'pert_dose':{'$gt':5}},{'sig_id':True},limit=1)
        #no dose min
        pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell},{'sig_id':True},limit=1)
        if pert_query:
            sigs.append(pert_query[0])
        # else:
        #     print 'no sig for ' + pert + ' in ' + cell
    with open(file1, 'a') as f:
        [f.write(x + '\n') for x in sigs]
metric = 'wtcs'
queryDir = work_dir + '/ctd2_sig_query'
if not os.path.exists(queryDir):
    os.mkdir(queryDir)
cmd = ' '.join(['rum -q local -f sig_query_tool',
         '--sig_id ' + file1,
         '--metric ' + metric,
         '--column_space full',
         '--out ' + queryDir,
         '--mkdir false',
         '--save_tail false'])
         # '--row_space bing', 
os.system(cmd)


metric = 'wtcs'
queryDir = work_dir + '/test_sig_query'
if not os.path.exists(queryDir):
    os.mkdir(queryDir)
file1 = '/xchip/cogs/hogstrom/analysis/summly/cp_class/sigs_combo.grp'
cmd = ' '.join(['rum -q local -f sig_query_tool',
         '--sig_id ' + file1,
         '--metric ' + metric,
         '--column_space full',
         '--out ' + queryDir,
         '--mkdir false',
         '--save_tail false'])

### run sig_summly with matlab 

# read in summly 
sigFrame = pd.read_csv(file1,header=None,names=['sigs'])
sigSer = sigFrame['sigs']
sigSet = set(sigSer)
sigSmSer = pd.Series(list(sigSet))
file2 = work_dir + '/' + 'ctd2_sig_set.grp'
sigSmSer.to_csv(file2,index=False)


map1 = { 'a': 1, 'b':2 }
# inv_map = {[v:k for k, v in map1.items()]}
inv_map = dict(zip(map1.values(), map1.keys()))




########################################
### read + explore summly results ######
########################################

# subset of molecular targets + cps
# bigGroups = ['AKT','BRAF','EGFR','HDAC','MAPK','MTOR','HSP90','PI3K','PPARG','TP53']
# brdGrpList = []
# grpToCp = {}
# for grp in bigGroups:
#     grpPerts = drugLabels['pert_id'][drugLabels['gene_dw_annot'] == grp]
#     grpToCp[grp] = list(grpPerts.values)
#     brdGrpList.extend(grpPerts.values)
# get list of cps in summly dir
basePath = work_dir + '/sig_query'
dateID = 'jul30/my_analysis.sig_summly_tool.2013073021435191'
summDir = '/'.join([basePath,dateID])
cpDirs = [f for f in os.listdir(summDir) if os.path.isdir(summDir+'/'+f)]

graphDir = work_dir + '/graph_out_4Sept'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
### examine each functional gorup
for grpGene in grpToCp:
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
    ## scratch
    # plt.tight_layout(pad=2, w_pad=4, h_pad=1.0)
    # plt.tight_layout()
    # plt.savefig(outF, bbox_inches='tight')
    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    # plt.savefig()
    
    # plt.savefig(outF, bbox_extra_artists=(lgd,), bbox_inches='tight')
    # fig.savefig(outF, bbox_inches='tight')


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
# basePath = work_dir + '/sig_query'
basePath = '/xchip/cogs/hogstrom/analysis/summly/cp_class/ctd2_sig_query'
dateID = 'aug01/my_analysis.sig_summly_tool.2013080119394091'
summDir = '/'.join([basePath,dateID])
cpDirs = [f for f in os.listdir(summDir) if os.path.isdir(summDir+'/'+f)]

graphDir = work_dir + '/ctd2_graph_out_4Sept'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)

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

### examine each functional gorup
sumScoreDict = {} #matrix of connections for each group
percSummlyDict = {}
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
    #take averages
    av_grp_sum_score = av_mtrx(grp_sum_score)
    av_grp_PercSummly = av_mtrx(grp_PercSummly)
    av_grp_rank = av_mtrx(grp_rank)
    # store averages
    sumScoreDict[grpGene] = grp_sum_score
    percSummlyDict[grpGene] = grp_PercSummly
    ### print group heatmap
    fig = plt.figure(1, figsize=(20, 8))
    plt.suptitle(grpGene + ' compound group',fontsize=14, fontweight='bold')
    plt.subplot(121)
    plt.title('percent summly rank')
    plt.imshow(av_grp_PercSummly,
            interpolation='nearest',
            cmap=matplotlib.cm.Greens_r,
            vmin=0, 
            vmax=1)
    ytcks = [pDescDict[x] for x in grp]
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
    outF = os.path.join(graphDir,grpGene + '_compound_group_heatmap.png')
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

