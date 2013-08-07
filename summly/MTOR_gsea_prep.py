'prep MTOR drugs for GSEA analysis'
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

gseaDir = work_dir + '/MTOR_cp_GSEA_lists'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
### examine each functional gorup
grpGene = 'MTOR'
grp = grpToCp[grpGene]
# matrices for grp connections
nGrp = len(grp)
grp_sum_score = np.zeros((nGrp,nGrp))
grp_PercSummly = np.zeros((nGrp,nGrp))
grp_rank = np.zeros((nGrp,nGrp))
for ibrd,brd in enumerate(grp):
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
    #write list of kd/oe genes to a file for gsea preranked
    outF = gseaDir + '/' + brd + '_KD_summly_result.txt'
    scoreSeries = pd.Series(cgsRes['sum_score'].values, index=cgsRes['pert_iname'])
    scoreSeries.to_csv(path=outF,sep='\t')



