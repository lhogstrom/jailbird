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

filePCLgrps = '/xchip/cogs/projects/pharm_class/pcl_shared_target.txt'
pclFrm = pd.io.parsers.read_csv(filePCLgrps,sep='\t')
#replace ugly characters
pclFrm['class'] = pclFrm['class'].str.replace("/","-")
pclFrm['class'] = pclFrm['class'].str.replace(" ","_")
pclFrm['class'] = pclFrm['class'].str.replace("&","_")
pclFrm['class'] = pclFrm['class'].str.replace("?","_")
pclFrm['class'] = pclFrm['class'].str.replace("(","_")
pclFrm['class'] = pclFrm['class'].str.replace(")","_")
pclFrm['class'] = pclFrm['class'].str.replace("'","_")

# TestGrps = ['ATP1A1','GSK3B','KDR']
grpDict = {}
for grp in pclFrm['class']:
    grpFrm = pclFrm[pclFrm['class'] == grp]
    grpBrds = grpFrm['pert_id'].values
    grpDict[grp] = list(grpBrds)


### class analysis
wkdir = '/xchip/cogs/projects/pharm_class/summly_Lass_Dec14'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
reload(pcla)
metric = 'wtcs'
po = pcla.PCLA(grpDict,
                    metric,
                    wkdir,
                    pairwise_prefix='pairwise_matrices',
                    rankpt_metric='mean_rankpt_4',
                    sum_score_metric='sum_score_4',
                    row_space = 'lm',
                    cell_match_mode=True)
po.get_inames()
po.load_summly_mtrx(skip_brd_duplicates=False)
po.test_groups(make_heatmaps=True,
            group_size_min=3)
po.make_summary_boxplot()
po.test_class_interrelatedness(make_heatmaps=True,
        make_boxplots=True,
        make_group_by_cp_mtrx=True)
po.inter_group_cluster(po.interGrpFrm)
po.inter_group_line_graph()

#don't want the same connections drive gorups over and over again.
#1) mask allBRD rnkpt matrix to threshold (rnkpt > 90). 
#2) make upper matrix
#3) if the same pert-to-pert connection, ignore
#2) across all PCLs, how many times is a specific drug-drug connection seen as driving a PCL



