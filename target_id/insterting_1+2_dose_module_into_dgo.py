import glob, HTML
import matplotlib.pyplot as plt
import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection as fdr
import cmap.util.mongo_utils as mu
from cmap.tools import sig_slice_tool
from cmap.io import gct,plategrp,rnk
import cmap.util.progress as progress
import subprocess
import cmap.util.tool_ops as to
import cmap.analytics.dgo as dgo
import cmap.analytics.dgo_oracle as dgo_oracle

import unittest
import numpy as np
from cmap.analytics import oracle
from cmap.io import queryresult

work_dir = '/xchip/cogs/projects/target_id/oracle_dose_5Aug'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)
#test query
pert_list = ['BRD-K69840642', 'BRD-K41859756']
pDescDict = {'BRD-K69840642':'ISOX','BRD-K41859756':'NVP-AUY922'}
targetDict = {'BRD-K69840642':['HSPA9'],'BRD-K41859756':['HSPA5']}

### dgo_oracle setup
dg = dgo_oracle.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# dg.get_sig_ids(genomic_pert='KD',targetDict_loaded=False,pert_list=pert_list,is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
dg.make_qres_obj(metric='spearman',gp_type='KD')
dg.test_known_connections(gp_type='KD',metric='spearman',
                            outName='apriori_connection_graphs',niter=10000,
                            conn_thresh=.01,connection_direction='both',
                            make_graphs=True,dose_graph=False)




### oracle scratch
direction = 'both'
q = targetDict.keys()[0]
grp1 = targetDict[q][0]
grpInnerMtch = [x for x in dg.ocl.col_groups.keys() if x.split('_')[0] == grp1]
g = grpInnerMtch[0]
rnk_vals = dg.ocl.query.pctrank.ix[dg.ocl.row_groups[q],dg.ocl.col_groups[g]]
score_vals = dg.ocl.query.score.ix[dg.ocl.row_groups[q],dg.ocl.col_groups[g]]
#calculate p-value
entry = dg.ocl.method_(dg.ocl.query,rnk_vals,q,g,direction)
# summary_list.append(entry)
rnk_vals = rnk_vals.unstack()
score_vals = score_vals.unstack()
rnk_vals = rnk_vals[rnk_vals.notnull()]
score_vals = score_vals[score_vals.notnull()]



        for q in targetDict:
            for grp1 in targetDict[q]:
                #find the indices matching first grouping term
                grpInnerMtch = [x for x in self.col_groups.keys() if x.split('_')[0] == grp1]
                for g in grpInnerMtch:
                    print "Connecting to group %s"  % str(g)
                    rnk_vals = self.query.pctrank.ix[self.row_groups[q],self.col_groups[g]]
                    score_vals = self.query.score.ix[self.row_groups[q],self.col_groups[g]]
                    #calculate p-value
                    entry = self.method_(self.query,rnk_vals,q,g,direction)
                    summary_list.append(entry)
                    rnk_vals = rnk_vals.unstack()
                    score_vals = score_vals.unstack()
                    rnk_vals = rnk_vals[rnk_vals.notnull()]
                    score_vals = score_vals[score_vals.notnull()]