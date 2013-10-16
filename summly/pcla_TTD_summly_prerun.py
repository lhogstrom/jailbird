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
import cmap.util.mongo_utils as mu

### Aravind's pcl class file ~300 cps
# drugFile = '/xchip/cogs/projects/pharm_class/pcl_classes.txt'
# drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
# #repalce ugly characters
# drugLabels['class'] = drugLabels['class'].str.replace("/","-")
# drugLabels['class'] = drugLabels['class'].str.replace(" ","_")
# drugLabels['class'] = drugLabels['class'].str.replace("&","_")
# drugLabels['class'] = drugLabels['class'].str.replace("?","_")
# drugLabels['class'] = drugLabels['class'].str.replace("(","_")
# drugLabels['class'] = drugLabels['class'].str.replace(")","_")
# # set up dictionary of compound classes
# grpSet = set(drugLabels['class'])
# grpToCp = {}
# for grp in grpSet:
#     grpPerts = drugLabels['pert_id'][drugLabels['class'] == grp]
#     grpToCp[grp] = list(grpPerts.values)
# # compound to group dict
# cpToGrp = {}
# for ibrd, brd in enumerate(drugLabels['pert_id']):
#     cpToGrp[brd] = drugLabels['class'][ibrd]
# inameDict = {}
# for ibrd,brd in enumerate(drugLabels['pert_id']):
#     inameDict[brd] = drugLabels.ix[ibrd]['pert_iname']

# ## TTD pcl class list - updated grpSetept/2013
drugFile = '/xchip/cogs/hogstrom/notes/TTD_annotations/TTD_targetID_pert_id_mapping_20130927.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
# replace ugly characters
drugLabels['Name'] = drugLabels['Name'].str.replace("/","-")
drugLabels['Name'] = drugLabels['Name'].str.replace(" ","_")
drugLabels['Name'] = drugLabels['Name'].str.replace("&","_")
drugLabels['Name'] = drugLabels['Name'].str.replace("?","_")
drugLabels['Name'] = drugLabels['Name'].str.replace("(","_")
drugLabels['Name'] = drugLabels['Name'].str.replace(")","_")
drugLabels['Name'] = drugLabels['Name'].str.replace("'","_")
# specify directionality
brdLstLst = drugLabels.pert_ids_merged.str.split('|').tolist()
drugLabels['pert_ids_merged'] = brdLstLst
grouped = drugLabels.groupby(['TTD_target_ID','Category'])
ttd_cp_dict = {}
for ig, group in enumerate(grouped):
    # if ig > 20: # shorten diction 
    #     continue
    groupName = group[1]['Name'].values[0]
    groupCat = group[1]['Category'].values[0]
    cpLstLst = group[1]['pert_ids_merged'].values
    cpLst = [item for sublist in cpLstLst for item in sublist] 
    cpLst = list(set(cpLst)) # make sure list is unique
    ttd_cp_dict[groupName+'-'+groupCat] = cpLst
#set inames
inameDict = {}
for x in drugLabels.iterrows():
    brds = x[1]['pert_ids_merged']
    iname = x[1]['Drug_name']
    if len(brds) > 1:
        for brd in brds:
            inameDict[brd] = iname
    else:
        inameDict[brds[0]] = iname


### drugBank list
drugFile = '/xchip/cogs/hogstrom/notes/TTD_annotations/drugbanktarget_pert_id.csv'
drugLabels = pd.io.parsers.read_csv(drugFile)
drugLabels['action_calc'] = drugLabels['action_calc'].str.replace("other/unknown","other")
# specify directionality
brdLstLst = drugLabels.pert_ids.str.split('|').tolist()
drugLabels['pert_ids'] = brdLstLst
# group by gene and action
grouped = drugLabels.groupby(['gene','action_calc'])
ttd_cp_dict = {}
for group in grouped:
    groupName = group[1]['gene'].values[0]
    groupCat = group[1]['action_calc'].values[0]
    cpLstLst = group[1]['pert_ids'].values
    cpLst = [item for sublist in cpLstLst for item in sublist] 
    cpLst = list(set(cpLst)) # make sure list is unique
    if groupCat == '-666':
        ttd_cp_dict[groupName] = cpLst
    else:
        ttd_cp_dict[groupName+'-'+groupCat] = cpLst
# group only be gene
grouped = drugLabels.groupby(['gene'])
gene_cp_dict = {}
for ig, group in enumerate(grouped):
    # if ig > 40:
    #     continue
    if group[0] == '-666':
        continue    
    groupName = group[1]['gene'].values[0]
    cpLstLst = group[1]['pert_ids'].values
    cpLst = [item for sublist in cpLstLst for item in sublist] 
    cpLst = list(set(cpLst)) # make sure list is unique
    gene_cp_dict[groupName] = cpLst


### class analysis
wkdir = '/xchip/cogs/projects/pharm_class/TTd_Oct16'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
reload(pcla)
metric = 'wtcs'
po = pcla.PCLA(ttd_cp_dict,    
                    metric,
                    wkdir,
                    pairwise_prefix='pairwise_matrices',
                    rankpt_metric='mean_rankpt_4',
                    sum_score_metric='sum_score_4',
                    row_space = 'lm',
                    cell_match_mode=True)
po.get_inames()
po.load_summly_mtrx()
po.test_groups(make_heatmaps=True,
            group_size_min=5)
po.make_summary_boxplot()

### Drug-KD connection
wkdir = '/xchip/cogs/projects/pharm_class/drugBank_Oct16'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
reload(pcla)
metric = 'wtcs'
po = pcla.PCLA(gene_cp_dict,    
                    metric,
                    wkdir,
                    pairwise_prefix='pairwise_matrices',
                    rankpt_metric='mean_rankpt_4',
                    sum_score_metric='sum_score_4',
                    row_space = 'lm',
                    cell_match_mode=True)
po.get_inames()
po.load_summly_mtrx()
po.check_shRNA_connection() 


# po2 = po
po.sum_scoreMtrx = po2.sum_scoreMtrx
po.rankptInFrm = po2.rankptInFrm
po.rankptMtrx = po2.rankptMtrx
po.cp_percent_summlyMtrx = po2.cp_percent_summlyMtrx
po.pclResultDict = po2.pclResultDict
po.brdTpDict = po2.brdTpDict

## count number of compounds with a drugBank TTD_target_ID
subList = []
for x in po.pclResultDict:
    subList.append(po.pclResultDict[x])
cpList = [item for sublist in subList for item in sublist]
len(set(cpList))


