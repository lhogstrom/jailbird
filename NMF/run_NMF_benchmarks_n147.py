#! /usr/bin/env python

'''
run NMF_benchmarks module

LH 02/2014
'''
import os
import cmap.analytics.NMF_benchmarks as nmfb
import pandas as pd

# nComponents = 9
# dimDict = {'A375_c9_lm_epsilon':'n1057x978',
# 'A549_c9_lm_epsilon':'n1405x978', # 
# 'HA1E_c9_lm_epsilon':'n1345x978',
# 'HCC515_c9_lm_epsilon':'n1205x978',
# 'HEPG2_c9_lm_epsilon':'n783x978',
# 'HT29_c9_lm_epsilon':'n899x978',
# 'MCF7_c9_lm_epsilon':'n1472x978',
# 'PC3_c9_lm_epsilon':'n1327x978',
# 'VCAP_c9_lm_epsilon':'n1359x978'}

nComponents = 20
dimDict = {'A375_c20_lm_epsilon':'n1057x978',
'A549_c20_lm_epsilon':'n1405x978', # 
'HA1E_c20_lm_epsilon':'n1345x978',
'HCC515_c20_lm_epsilon':'n1205x978',
'HEPG2_c20_lm_epsilon':'n783x978',
'HT29_c20_lm_epsilon':'n899x978',
'MCF7_c20_lm_epsilon':'n1472x978',
'PC3_c20_lm_epsilon':'n1327x978',
'VCAP_c20_lm_epsilon':'n1359x978'}
# cellList = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines

###############################
### NMF benchmarks on drugs ###
###############################

wkdir = '/xchip/cogs/projects/NMF/cpd_groups_n147'
for prefix in dimDict:
    dim = dimDict[prefix]
    # if cell == 'A549':
    #     continue # skip failed job
    # this W matrix will be used
    path1 = wkdir + '/' + prefix
    prefix1 = 'clique_compound_classes_' + dim
    outdir = path1 + '/benchmark_graphs'
    source_dir = path1
    Hfile = prefix1 + '.H.k' + str(nComponents) + '.gct'
    # WFile = prefix1 + '.W.k' + str(nComponents) + '.gct'
    # MI_file_component = prefix1 + '.MI.k' + str(nComponents) + '.gct'
    # MI_file_inspace = prefix1 + '.MI.input_space.gct'
    anntFile = 'clique_compound_classes.v2.txt'
    # groupFile = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/n69_drug_targets.gmt'
    groupFile = '/xchip/cogs/projects/pharm_class/rnwork/cliques/cpd_groups_n147.gmt'
    # run NMF module 
    reload(nmfb)
    self = nmfb.NMFresult(source_dir)
    self.set_output_dir(out=outdir)
    self.load_NMF_H_matrix(Hfile)
    self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
    self.group_component_maps(match_field='pert_id')
    # MI - 9 components
    # self.load_MI_matrix(MI_file_component)
    # self.MI_pairwise_comp(self.mi,match_field='pert_id')
    # self.intra_group_boxplot(space_name='c9')
    # self.boxplot_with_null(space_name='c9')
    # # MI - LM
    # self.load_MI_matrix(MI_file_inspace)
    # self.MI_pairwise_comp(self.mi,match_field='pert_id')
    # self.intra_group_boxplot(space_name='LM')
    # self.boxplot_with_null(space_name='LM')
