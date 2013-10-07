'''
Contains code illustrating usage examples of the pcla class - pharmacological class analyzer
'''
import cmap.analytics.pcla as pcla
import numpy as np
import os
import cmap.io.gmt as gmt
import cmap
import pandas as pd
import matplotlib.pyplot as plt

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

inameDict = {}
for ibrd,brd in enumerate(drugLabels['pert_id']):
    inameDict[brd] = drugLabels.ix[ibrd]['pert_iname']

### cell line match mode
wkdir = '/xchip/cogs/sig_tools/sig_summly/pcl/dw_sep_24'
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
po.get_sig_ids()
# po.run_summly(rerun_mode=False)
# summPath = po.out + '/summly_out/sep11'
summPath = '/xchip/cogs/projects/connectivity/summly/matched/src'
po.make_summly_path_dict(summPath)
# po.run_summly(rerun_mode=True)
# po.make_summly_path_dict(summPath_nMtch)
po.get_inames()
# po.test_groups(make_heatmaps=True,
#         group_size_min=3,
#         sum_score_metric='sum_score_4',
#         rankpt_metric='mean_rankpt_4')
# po.cluster_all_cps()
# po.make_summary_boxplot()
# po.cluster_all_cps(make_heatmaps=True,
#         sum_score_metric='sum_score_4',
#         rankpt_metric='mean_rankpt_4')
po.test_class_interrelatedness(make_heatmaps=False,
                            make_boxplots=False, 
                            rankpt_metric='mean_rankpt_4')
# po.inter_group_cluster(self.interGrpFrm)

