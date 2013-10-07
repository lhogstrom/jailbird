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
inameDictDOS = {}
for brd in dosbioSet:
    if brd in sumFrm.index:
        grpToCp['DOS'].append(brd)
        inameDictDOS[brd] = brd

### cell line match mode
wkdir = '/xchip/cogs/sig_tools/sig_summly/dosbio/Match_analysis_Sept16'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
reload(pcla)
metric = 'wtcs'
poDOS = pcla.PCLA(grpToCp,
                    metric,
                    wkdir,
                    summly_out_prefix='summly_out',
                    pairwise_prefix='pairwise_matrices',
                    cell_match_mode=True, 
                    row_space = 'lm')
poDOS.get_sig_ids()
# po.run_summly(rerun_mode=False)
# summPath = po.out + '/summly_out/sep11'
summPath = '/xchip/cogs/sig_tools/sig_summly/dosbio/summly_out/sep06'
poDOS.make_summly_path_dict(summPath)
# # po.run_summly(rerun_mode=True)
# # po.make_summly_path_dict(summPath_nMtch)
# po.get_inames()
# po.test_groups(make_heatmaps=True,
#         group_size_min=3,
#         sum_score_metric='sum_score_4',
#         rankpt_metric='mean_rankpt_4')
# po.cluster_all_cps()
# po.make_summary_boxplot()

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
grpToCp = {}
for grp in grpSet:
    grpPerts = drugLabels['pert_id'][drugLabels['class'] == grp]
    grpToCp[grp] = list(grpPerts.values)
# compound to group dict
cpToGrp = {}
for ibrd, brd in enumerate(drugLabels['pert_id']):
    cpToGrp[brd] = drugLabels['class'][ibrd]

inameDictPCL = {}
for ibrd,brd in enumerate(drugLabels['pert_id']):
    inameDictPCL[brd] = drugLabels.ix[ibrd]['pert_iname']


### cell line match mode
wkdir = '/xchip/cogs/sig_tools/sig_summly/pcl/Match_Sept16'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
reload(pcla)
metric = 'wtcs'
poPCL = pcla.PCLA(grpToCp,
                    metric,
                    wkdir,
                    summly_out_prefix='summly_out',
                    pairwise_prefix='pairwise_matrices',
                    cell_match_mode=True, 
                    row_space = 'lm')
poPCL.get_sig_ids()
# po.run_summly(rerun_mode=False)
# summPath = po.out + '/summly_out/sep11'
summPath = '/xchip/cogs/sig_tools/sig_summly/pcl/summly_out/sep06'
poPCL.make_summly_path_dict(summPath)

### merge dictionaries of PCL and DOS analysis 
pathDictAll = dict(poDOS.cpPathDict.items() + poPCL.cpPathDict.items())
grpToCpAll = dict(poDOS.pclDict.items() + poPCL.pclDict.items())
cpSetAll = poDOS.cpSet.union(poPCL.cpSet)
inameDictAll = dict(inameDictDOS.items() + inameDictPCL.items())

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
                    pairwise_prefix='pairwise_matrices',
                    cell_match_mode=True, 
                    row_space = 'lm')
# po.get_sig_ids()
po.cpPathDict = pathDictAll
po.cpSet = cpSetAll
po.pclDict = grpToCpAll
po.inameDict = inameDictAll
# po.run_summly(rerun_mode=False)
# summPath = po.out + '/summly_out/sep11'
# summPath = '/xchip/cogs/sig_tools/sig_summly/dosbio/summly_out/sep06'
# po.make_summly_path_dict(summPath)
# po.run_summly(rerun_mode=True)
# po.make_summly_path_dict(summPath_nMtch)
# po.get_inames()
po.test_groups(make_heatmaps=False,
        group_size_min=3,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4')
# po.cluster_all_cps()
po.make_summary_boxplot()




### scratch on DOS summly anlaysis
self = po
make_heatmaps=True
group_size_min=3
sum_score_metric='sum_score_4'
rankpt_metric='mean_rankpt_4'


dos_rnkpt = {}
for doscp in self.pclDict['DOS']
    dos_rnkpt_dist = {}
    for grpName in self.pclDict: # top connecting cps?
        if (grpName == '-666') or (grpName == '-666'):
            continue
        grp = self.pclDict[grpName]
        grp.append(doscp)
        if not grp: # skip if grp is empty
            continue
        grp = [cp for cp in grp if cp in self.cpPathDict] # leave out compounds that don't have summly data
        if len(grp) < group_size_min:
            continue
        grpInames = [self.inameDict[cp] for cp in grp] # inames
        grpZip = zip(*[grp,grpInames])
        # matrices for group connections
        if self.cell_match_mode:
            [grp_sum_score, 
            grp_rankpt, 
            grp_PercSummly] = self.pairwise_calc_match(grp,
                                            sum_score_metric,
                                            rankpt_metric)
        else:
            [grp_sum_score, 
            grp_rankpt, 
            grp_PercSummly] = self.pairwise_calc_nonmatch(grp,
                                            sum_score_metric,
                                            rankpt_metric)
        ### take averages of the upper and lower matrix segments
        av_grp_sum_score = self.av_mtrx(grp_sum_score)
        av_grp_rankpt = self.av_mtrx(grp_rankpt)
        av_grp_PercSummly = self.av_mtrx(grp_PercSummly)
        dos_rnkpt_dist[grpName] = av_grp_rankpt[:,-1]
    dos_rnkpt[doscp] = dos_rnkpt_dist
self.dos_rnkpt_by_goup = dos_rnkpt


