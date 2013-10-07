'''
Contains code illustrating usage examples of the pcla class - pharmacological class analyzer
'''
import cmap.analytics.pcla as pcla
import cmap.util.mongo_utils as mu
import numpy as np
import os
import cmap.io.gmt as gmt
import cmap
import pandas as pd
import matplotlib.pyplot as plt

### get DOS BIO cps from mongo
CM = mu.CMapMongo()
pert_List = CM.find({'sig_id':{'$regex':'DOSBIO'},'pert_iname':{'$regex':'BRD'}},{'pert_id':True})
dosbioSet = set(pert_List)

# load in sumSpace file
sumSpaceFile = '/xchip/cogs/sig_tools/sig_query/results/a2_internal_wtcs.lm/query_info_n73597.txt'
summSpace = pd.read_csv(sumSpaceFile,sep='\t')

#make smaller summspace Frame
sumFrm = summSpace.copy()
sumFrm.index = summSpace['pert_id']
# identify DOSBIOs in summ space
grpToCp = {}
grpToCp['DOS'] = []
inameDict = {}
for brd in dosbioSet:
    if brd in sumFrm.index:
        grpToCp['DOS'].append(brd)
        inameDict[brd] = brd

### pcla for knwon groups
drugFile = '/xchip/cogs/sig_tools/sig_summly/pcl/pcl_classes.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
#repalce ugly characters
drugLabels['class'] = drugLabels['class'].str.replace("/","-")
drugLabels['class'] = drugLabels['class'].str.replace(" ","_")
drugLabels['class'] = drugLabels['class'].str.replace("&","_")
drugLabels['class'] = drugLabels['class'].str.replace("?","_")
drugLabels['class'] = drugLabels['class'].str.replace("(","_")
drugLabels['class'] = drugLabels['class'].str.replace(")","_")

# set up dictionary of compound classes
grpSet = set(drugLabels['class'])
# grpToCp = {}
for grp in grpSet:
    grpPerts = drugLabels['pert_id'][drugLabels['class'] == grp]
    grpToCp[grp] = list(grpPerts.values)

# inameDict = {}
for ibrd,brd in enumerate(drugLabels['pert_id']):
    inameDict[brd] = drugLabels.ix[ibrd]['pert_iname']

### merge dictionaries of PCL and DOS analysis 
# pathDictAll = dict(poDOS.cpPathDict.items() + poPCL.cpPathDict.items())
# grpToCpAll = dict(poDOS.pclDict.items() + poPCL.pclDict.items())
# cpSetAll = poDOS.cpSet.union(poPCL.cpSet)
# inameDictAll = dict(inameDictDOS.items() + inameDictPCL.items())
# po.cpPathDict = pathDictAll
# po.cpSet = cpSetAll
# po.pclDict = grpToCpAll
# po.inameDict = inameDictAll


### combined PCLA 
wkdir = '/xchip/cogs/sig_tools/sig_summly/dosbio/Match_analysis_Sept16'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
reload(pcla)
metric = 'wtcs'
po = pcla.PCLA(grpToCp,
                    metric,
                    wkdir,
                    summly_out_prefix='summly_out',
                    pairwise_prefix='pairwise_matrices2',
                    cell_match_mode=True, 
                    row_space = 'lm')
po.get_sig_ids(write_grps=False)
# po.run_summly(rerun_mode=False)
# summPath = po.out + '/summly_out/sep11'
summPath = '/xchip/cogs/projects/connectivity/summly/matched/src'
po.make_summly_path_dict(summPath)
# po.run_summly(rerun_mode=True)
# po.make_summly_path_dict(summPath_nMtch)
po.get_inames()
po.test_DOS_queries(make_heatmaps=True,
        make_boxplots=True, 
        rnkpt_thresh=90,       
        group_size_min=3,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4')
# po.test_groups(make_heatmaps=False,
#         group_size_min=3,
#         sum_score_metric='sum_score_4',
#         rankpt_metric='mean_rankpt_4')
# po.make_summary_boxplot()
# po.cluster_all_cps()

# self = po
# make_heatmaps=True
# group_size_min=3
# sum_score_metric='sum_score_4'
# rankpt_metric='mean_rankpt_4'
