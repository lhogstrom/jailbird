'''
Contains code illustrating usage examples of the pcla class - pharmacological class analyzer
'''
import cmap.analytics.pcla as pcla
import numpy as np
import os
import cmap.io.gmt as gmt
import cmap
import pandas as pd

summPath = '/xchip/cogs/sig_tools/sig_summly/pcl/summly_out/sep06'
wkdir = '/xchip/cogs/sig_tools/sig_summly/pcl/graphs_Sept11'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

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

reload(pcla)
gp_type = 'KD'
metric = 'wtcs'
pclaObj = pcla.PCLA(grpToCp,  
                    metric, 
                    wkdir)
pclaObj.get_sig_ids()
pclaObj.run_summly(cell_match_mode='true')
pclaObj.make_summly_path_dict(summPath)
pclaObj.inameDict = inameDict #make this part of the tool
pclaObj.test_groups(make_heatmaps=False,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4')
pclaObj.make_summary_boxplot()

reload(pcla)
pclaObj2 = pcla.PCLA(grpToCp,  
                    metric, 
                    wkdir)
pclaObj2.sumScoreDict = pclaObj.sumScoreDict
pclaObj2.rnkptDict = pclaObj.rnkptDict
pclaObj2.percSummlyDict = pclaObj.percSummlyDict
pclaObj2.make_summary_boxplot()


