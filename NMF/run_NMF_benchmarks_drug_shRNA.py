#! /usr/bin/env python

'''
run NMF_benchmarks module

LH 02/2014
'''
import os
import cmap.analytics.NMF_benchmarks as nmfb
import pandas as pd

### directory information
nComponents = 9
shRNA_dimDict = {'A375_shRNA_c9_lm_epsilon':'n279x978',
'A549_shRNA_c9_lm_epsilon':'n258x978', # 
'HA1E_shRNA_c9_lm_epsilon':'n267x978',
'HCC515_shRNA_c9_lm_epsilon':'n232x978',
'HEPG2_shRNA_c9_lm_epsilon':'n252x978',
'HT29_shRNA_c9_lm_epsilon':'n259x978',
'MCF7_shRNA_c9_lm_epsilon':'n400x978',
'PC3_shRNA_c9_lm_epsilon':'n441x978',
'VCAP_shRNA_c9_lm_epsilon':'n448x978'}
cp_dimDict = {'A375_drug_c9_lm_epsilon':'n340x978',
'A549_drug_c9_lm_epsilon':'n345x978', # 
'HA1E_drug_c9_lm_epsilon':'n380x978',
'HCC515_drug_c9_lm_epsilon':'n364x978',
'HEPG2_drug_c9_lm_epsilon':'n268x978',
'HT29_drug_c9_lm_epsilon':'n301x978',
'MCF7_drug_c9_lm_epsilon':'n450x978',
'PC3_drug_c9_lm_epsilon':'n428x978',
'VCAP_drug_c9_lm_epsilon':'n380x978'}
# put into dataframe
shRNASer = pd.Series(shRNA_dimDict)
cpSer = pd.Series(cp_dimDict)
cellList = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
# make data frames of the meta data
cpFrm = pd.DataFrame(data={'cp_dim':cpSer.values,
                    'cp_name':cpSer.index,
                    'shRNA_dim':shRNASer.values,
                    'shRNA_name':shRNASer.index},
                    index=cellList)

#############################
### shRNA of drug targets ###
#############################

wkdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA'
for cell in cpFrm.index:
    # this W matrix will be used
    path1 = wkdir + '/' + cpFrm.ix[cell,'shRNA_name']
    prefix1 = 'shRNA_drug_target_genes_' + cpFrm.ix[cell,'shRNA_dim']   
    # new H matrix will match the above W matrix
    # path2 = wkdir + '/' + cpFrm.ix[cell,'cp_name']
    # prefix2 = 'clique_compound_classes_' + cpFrm.ix[cell,'cp_dim']
    # nFolder = 'shRNA_projections'
    # nFPath = path2 + '/' +nFolder
    outdir = path1 + '/benchmark_graphs'
    source_dir = path1
    Hfile = prefix1 + '.H.k' + str(nComponents) + '.gct'
    WFile = prefix1 + '.W.k' + str(nComponents) + '.gct'
    MI_file_component = prefix1 + '.MI.k' + str(nComponents) + '.gct'
    MI_file_inspace = prefix1 + '.MI.input_space.gct'
    anntFile = 'shRNA_drug_target_genes.v2.txt'
    groupFile = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/n69_shRNA_sig_ids.gmt'
    # run NMF module 
    reload(nmfb)
    self = nmfb.NMFresult(source_dir)
    self.set_output_dir(out=outdir)
    self.load_NMF_H_matrix(Hfile)
    self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
    self.group_component_maps(match_field='pert_iname')
    # MI - 9 components
    # self.load_MI_matrix(MI_file_component)
    # self.MI_pairwise_comp(self.mi,match_field='pert_iname')
    # self.intra_group_boxplot(space_name='c9')
    # self.boxplot_with_null(space_name='c9')
    # # MI - LM
    # self.load_MI_matrix(MI_file_inspace)
    # self.MI_pairwise_comp(self.mi,match_field='pert_iname')
    # self.intra_group_boxplot(space_name='LM')
    # self.boxplot_with_null(space_name='LM')

####################
### NMF on drugs ###
####################

wkdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA'
for cell in cpFrm.index:
    # if cell == 'A549':
    #     continue # skip failed job
    # this W matrix will be used
    path1 = wkdir + '/' + cpFrm.ix[cell,'cp_name']
    prefix1 = 'clique_compound_classes_' + cpFrm.ix[cell,'cp_dim']
    outdir = path1 + '/benchmark_graphs'
    source_dir = path1
    Hfile = prefix1 + '.H.k' + str(nComponents) + '.gct'
    WFile = prefix1 + '.W.k' + str(nComponents) + '.gct'
    MI_file_component = prefix1 + '.MI.k' + str(nComponents) + '.gct'
    MI_file_inspace = prefix1 + '.MI.input_space.gct'
    anntFile = 'clique_compound_classes.v2.txt'
    groupFile = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/n69_drug_targets.gmt'
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

################################################ 
### drug signatures projected in shRNA space ###
################################################ 

wkdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA'
for cell in cpFrm.index:
    if cell == 'A549':
        continue # skip failed job
    # this W matrix will be used
    path1 = wkdir + '/' + cpFrm.ix[cell,'shRNA_name'] + '/drug_projections'
    prefix1 = 'shRNA_drug_target_genes_' + cpFrm.ix[cell,'shRNA_dim']   
    # new H matrix will match the above W matrix
    path2 = wkdir + '/' + cpFrm.ix[cell,'cp_name']
    prefix2 = 'clique_compound_classes_' + cpFrm.ix[cell,'cp_dim']
    nFolder = 'shRNA_projections'
    nFPath = path2 + '/' +nFolder
    outdir = path1 + '/benchmark_graphs'
    source_dir = path1
    # Hfile = prefix +'_' + dim + '.H.k' + str(nComponents) + '.gct'
    Hfile = prefix1 + '.drug_projections.H.k' + str(nComponents) + '.gct'
    anntFile = 'shRNA_drug_target_genes.v2.txt'
    # groupFile = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/n69_drug_targets.gmt'
    groupFile = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/n69_shRNA_sig_ids.gmt'
    # run NMF module 
    reload(nmfb)
    self = nmfb.NMFresult(source_dir)
    self.set_output_dir(out=outdir)
    self.load_NMF_H_matrix(Hfile)
    self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
    self.group_component_maps(match_field='pert_iname')


################################################ 
### shRNA signatures projected in drug space ###
################################################ 

# output directory
outdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/VCAP_shRNA_c9_lm_epsilon/drug_projections/benchmark_graphs'
### Define Inputs
source_dir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/VCAP_shRNA_c9_lm_epsilon' # NMF result directory
nComponents = 9
prefix = 'VCAP_shRNA_c9_lm_epsilon'
dim = 'n448x978'
# Hfile = 'shRNA_drug_target_genes_' + dim + '.H.k' + str(nComponents) + '.gct'
Hfile = 'shRNA_drug_target_genes_' + dim + '.drug_projections.H.k' + str(nComponents) + '.gct'
WFile = 'shRNA_drug_target_genes_' + dim + '.W.k' + str(nComponents) + '.gct'
MI_file_component = 'shRNA_drug_target_genes_n448x978.MI.k' + str(nComponents) + '.gct'
MI_file_inspace = 'shRNA_drug_target_genes_n448x978.MI.input_space.gct'
anntFile = 'shRNA_drug_target_genes.v2.txt'
groupFile = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/n69_shRNA_sig_ids.gmt'

# run NMF module 
reload(nmfb)
self = nmfb.NMFresult(source_dir)
self.set_output_dir(out=outdir)
self.load_NMF_H_matrix(Hfile)
self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
self.group_component_maps(match_field='pert_iname')
# MI - 9 components
self.load_MI_matrix(MI_file_component)
self.MI_pairwise_comp(self.mi,match_field='pert_iname')
self.intra_group_boxplot(space_name='c9')
self.boxplot_with_null(space_name='c9')
# MI - LM
self.load_MI_matrix(MI_file_inspace)
self.MI_pairwise_comp(self.mi,match_field='pert_iname')
self.intra_group_boxplot(space_name='LM')
self.boxplot_with_null(space_name='LM')
