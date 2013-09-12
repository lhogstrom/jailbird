'''
Contains code illustrating usage examples of the pcla class - pharmacological class analyzer
'''
import cmap.analytics.pcla as pcla
import numpy as np
import os
import cmap.io.gmt as gmt
import cmap
import pandas as pd

wkdir = '/xchip/cogs/sig_tools/sig_summly/pcl/TTD_Sept11'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

drugFile = '/xchip/cogs/projects/cp_annot/ttd_cmap_export_filtered_20130910.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')

brdLstLst = drugLabels.pert_ids_merged.str.split('|').tolist()
drugLabels['pert_ids_merged'] = brdLstLst
grouped = drugLabels.groupby('TTD_target_ID')
grouped.groups

ttd_cp_dict = {}
for group in grouped:
    groupName = group[1]['Name'].values[0]
    cpLstLst = group[1]['pert_ids_merged'].values
    cpLst = [item for sublist in cpLstLst for item in sublist] 
    cpLst = list(set(cpLst))# make sure list is unique
    ttd_cp_dict[groupName] = cpLst

inameDict = {}
for x in drugLabels.iterrows():
    brds = x[1]['pert_ids_merged']
    iname = x[1]['Drug_name']
    for brd in brds:
        inameDict[brd] = iname

reload(pcla)
gp_type = 'KD'
metric = 'wtcs'
pclaObj = pcla.PCLA(ttd_cp_dict,  
                    metric, 
                    wkdir)
pclaObj.get_sig_ids()
pclaObj.run_summly(cell_match_mode='true')
summPath = pclaObj.out + '/sept11'
pclaObj.make_summly_path_dict(summPath)
pclaObj.inameDict = inameDict #make this part of the tool
pclaObj.test_groups(make_heatmaps=True,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4')
pclaObj.make_summary_boxplot()
