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

drugFile = '/xchip/cogs/projects/pharm_class/pcl_classes.txt'
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
wkdir = '/xchip/cogs/projects/pharm_class/Match_Sept27'
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
summPath = '/xchip/cogs/data/rnwork/batch_summly/summly_lm50'
po.make_rankpt_Mtrx(summPath)

po.get_sig_ids()
# po.run_summly(rerun_mode=False)
# summPath = po.out + '/summly_out/sep11'
# summPath = '/xchip/cogs/projects/connectivity/summly/matched/src'

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
po.cluster_all_cps(make_heatmaps=True,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4')
po.test_class_interrelatedness(make_heatmaps=True,
                            make_boxplots=True, 
                            make_group_by_cp_mtrx=True,
                            rankpt_metric='mean_rankpt_4')

### non match mode
wkdir = '/xchip/cogs/sig_tools/sig_summly/pcl/nonMatch_Sept14'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
reload(pcla)
metric = 'wtcs'
po = pcla.PCLA(grpToCp,
                    metric,
                    wkdir,
                    summly_out_prefix='summly_out',
                    pairwise_prefix='pairwise_matrices',
                    cell_match_mode=False, 
                    row_space = 'lm')
po.get_sig_ids()
# po.run_summly(rerun_mode=False)
# summPath = po.out + '/summly_out/sep11'
summPath_nMtch = '/xchip/cogs/sig_tools/sig_summly/pcl/summly_out_no_match/sep10'
po.make_summly_path_dict(summPath_nMtch)
po.run_summly(rerun_mode=True)
# po.make_summly_path_dict(summPath_nMtch)
po.get_inames()
po.test_groups(make_heatmaps=True,
        group_size_min=3,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4')
po.make_summary_boxplot()


#check which brds have fully completed brds
for cp in po.isSummSpace:
    sigs = po.sigIDdict[cp]
    for sig in sigs:
        jobIncomplete = False
        if not po.cpPathDict[cp].has_key(sig):
            print cp + ' : ' + sig + ' - incomplete '
            jobIncomplete = True
        else:
            path1 = cpPathDict[cp][sig]
            sumFile = path1 + '/'+sig + '_summly.txt'
            if not os.path.isfile(sumFile):
                jobIncomplete = True

### cp Path edit
cp = 'BRD-K08806317'
sig = 'CPC003_HA1E_6H_BRD-K08806317-050-03-6_10'
for cp in cpPathDict:
    sigs = self.sigIDdict[cp]
    for sig in sigs:
        jobIncomplete = False
        if not cpPathDict[cp].has_key(sig):
            jobIncomplete = True
        else:
            path1 = cpPathDict[cp][sig]
            sumFile = path1 + '/'+sig + '_summly.txt'
            if not os.path.isfile(sumFile):
                jobIncomplete = True
        if jobIncomplete:
            jobFalureList.append(cp)
            #remove direcotry for failed jobs
            if cpPathDict.has_key(cp):
                if cpPathDict[cp].has_key(sig):
                    cpDir = cpPathDict[cp][sig]
                    sDir = os.path.dirname(cpDir) #job direcotry
                    try:
                        print 'incomplete deleting: ' + sDir
                        shutil.rmtree(sDir)
                    except OSError:
                        continue

# list of operations to standardize:
# summly_handler tool = analogous to roger's querier?
# which sig_ids to use for a pert if there are more than one per cell line
# run summly
# manage paths to summly outputs 
# check that outputs have completed
# loading summly results
# combining the resulting _summly.txt files for non-matched mode

# would all of this be irreleant if the whole summly space pairwise matrix was calculated

self = po
allbrds = self.cpPathDict.keys()
grp = allbrds
rankptMtrx = po.whole_rankpt
percSummlyMtrx = po.whole_PercSummly

make_heatmaps=True
make_boxplots=True 
rnkpt_thresh=90     
group_size_min=3
sum_score_metric='sum_score_4'
rankpt_metric='mean_rankpt_4'

group_name = 'all_compounds'
grp = clust_order.values
rankpt_mtrx = rankptClustered.values
sumRank_mtrx = rankptClustered.values
out = outF


