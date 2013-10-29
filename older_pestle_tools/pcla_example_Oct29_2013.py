'''
Contains code illustrating usage examples of the pcla class - pharmacological class analyzer
'''
import cmap.analytics.pcla as pcla
import numpy as np
import os
import cmap.io.gmt as gmt
import cmap
import pandas as pd

#set working directory
wkdir = os.path.join(os.path.dirname(cmap.__file__), 'examples')
out = os.path.join(wkdir, 'pcla_example')

#define drug classes with members
pclDict = {'tubulin': ['BRD-K12539581',
                      'BRD-K35960502',
                      'BRD-A55594068',
                      'BRD-A76528577',
                      'BRD-K02407574',
                      'BRD-K11072542',
                      'BRD-K26997899',
                      'BRD-K32828673',
                      'BRD-K51318897',
                      'BRD-K77987382',
                      'BRD-K79131256'],
'microtubule-inhibitor': ['BRD-K12539581',
                      'BRD-A23723433',
                      'BRD-A60414806',
                      'BRD-A22783572',
                      'BRD-K08273968',
                      'BRD-K42125900']}
metric = 'wtcs'
### initiate the pcla object 
po = pcla.PCLA(pclDict,
            metric,
            out,
            summly_out_prefix='summly_out',
            pairwise_prefix='pairwise_matrices_by_Category',
            cell_match_mode=True, 
            row_space = 'lm')
po.get_sig_ids()
po.run_summly(rerun_mode=False)
summPath = po.out + '/summly_out/sep11'
po.make_summly_path_dict(summPath)
po.get_inames()
po.test_groups(make_heatmaps=False,
        group_size_min=15,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4')
po.make_summary_boxplot()
