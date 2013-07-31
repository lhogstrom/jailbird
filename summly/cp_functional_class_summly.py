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

allCtd2 = drugLabels['pert_id']
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
bigGroups = ['AKT','BRAF','EGFR','HDAC','MAPK','MTOR','HSP90','PI3K','PPARG','TP53']
brdGrpList = []
grpToCp = {}
for grp in bigGroups:
    grpPerts = drugLabels['pert_id'][drugLabels['gene_dw_annot'] == grp]
    grpToCp[grp] = list(grpPerts.values)
    brdGrpList.extend(grpPerts.values)
# get list of cps in summly dir
basePath = work_dir + '/sig_query'
dateID = 'jul30/my_analysis.sig_summly_tool.2013073021435191'
summDir = '/'.join([basePath,dateID])
cpDirs = [f for f in os.listdir(summDir) if os.path.isdir(summDir+'/'+f)]

graphDir = work_dir + '/graph_out'
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
    # print group heatmap
    plt.imshow(grp_rank,interpolation='nearest',cmap=matplotlib.cm.RdBu_r)
    plt.xticks(np.arange(len(grp)), grp,rotation=45)
    # ytcks = [pDescDict[x] for x in avicinsBrds]
    plt.yticks(np.arange(len(grp)),grp)
    plt.title(grpGene + ' group connections - sum_score')
    plt.colorbar()
    plt.savefig(os.path.join(graphDir,grpGene + '_compound_group_heatmap.png'))
    plt.close()
