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
drugLabels['Name'] = drugLabels['Name'].str.replace("/","-")
drugLabels['Name'] = drugLabels['Name'].str.replace(" ","_")
drugLabels['Name'] = drugLabels['Name'].str.replace("&","_")
drugLabels['Name'] = drugLabels['Name'].str.replace("?","_")
drugLabels['Name'] = drugLabels['Name'].str.replace("(","_")
drugLabels['Name'] = drugLabels['Name'].str.replace(")","_")


# all members of a family
# brdLstLst = drugLabels.pert_ids_merged.str.split('|').tolist()
# drugLabels['pert_ids_merged'] = brdLstLst
# grouped = drugLabels.groupby('TTD_target_ID')
# grouped.groups
# ttd_cp_dict = {}
# for group in grouped:
#     groupName = group[1]['Name'].values[0]
#     cpLstLst = group[1]['pert_ids_merged'].values
#     cpLst = [item for sublist in cpLstLst for item in sublist] 
#     cpLst = list(set(cpLst))# make sure list is unique
#     ttd_cp_dict[groupName] = cpLst

# specify directionality
brdLstLst = drugLabels.pert_ids_merged.str.split('|').tolist()
drugLabels['pert_ids_merged'] = brdLstLst
grouped = drugLabels.groupby(['TTD_target_ID','Category'])
ttd_cp_dict = {}
for group in grouped:
    groupName = group[1]['Name'].values[0]
    groupCat = group[1]['Category'].values[0]
    cpLstLst = group[1]['pert_ids_merged'].values
    cpLst = [item for sublist in cpLstLst for item in sublist] 
    cpLst = list(set(cpLst)) # make sure list is unique
    ttd_cp_dict[groupName+'-'+groupCat] = cpLst

inameDict = {}
for x in drugLabels.iterrows():
    brds = x[1]['pert_ids_merged']
    iname = x[1]['Drug_name']
    for brd in brds:
        inameDict[brd] = iname

reload(pcla)
gp_type = 'KD'
metric = 'wtcs'
po = pcla.PCLA(ttd_cp_dict,  
                    metric, 
                    wkdir,
                    summly_out_prefix='summly_out',
                    pairwise_prefix='pairwise_matrices_by_Category',
                    cell_match_mode=True, 
                    row_space = 'lm')
# po.get_sig_ids()
# po.run_summly(rerun_mode=False)
summPath = po.out + '/summly_out/sep11'
po.make_summly_path_dict(summPath)
# po.run_summly(rerun_mode=True)
po.make_summly_path_dict(summPath)
po.inameDict = inameDict #make this part of the tool
po.test_groups(make_heatmaps=False,
        group_size_min=15,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4')
po.make_summary_boxplot()


