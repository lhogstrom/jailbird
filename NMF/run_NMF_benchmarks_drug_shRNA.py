#! /usr/bin/env python

'''
run NMF_benchmarks module

LH 02/2014
'''
import os
import cmap.analytics.NMF_benchmarks as nmfb

#############################
### shRNA of drug targets ###
#############################

# output directory
outdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/VCAP_shRNA_c9_lm_epsilon/benchmark_graphs'
### Define Inputs
source_dir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/VCAP_shRNA_c9_lm_epsilon' # NMF result directory
nComponents = 9
prefix = 'VCAP_shRNA_c9_lm_epsilon'
dim = 'n448x978'
Hfile = 'shRNA_drug_target_genes_' + dim + '.H.k' + str(nComponents) + '.gct'
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
# self.group_component_maps(match_field='pert_iname')
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

####################
### NMF on drugs ###
####################

# output directory
outdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/VCAP_drug_c9_lm_epsilon/benchmark_graphs'
### Define Inputs
source_dir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/VCAP_drug_c9_lm_epsilon' # NMF result directory
nComponents = 9
prefix = 'clique_compound_classes'
dim = 'n380x978'
Hfile = prefix +'_' + dim + '.H.k' + str(nComponents) + '.gct'
WFile = prefix +'_' + dim + '.W.k' + str(nComponents) + '.gct'
MI_file_component = prefix + '_n448x978.MI.k' + str(nComponents) + '.gct'
MI_file_inspace = prefix + '_n448x978.MI.input_space.gct'
anntFile = prefix +'.v2.txt'
groupFile = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/n69_drug_targets.gmt'

# run NMF module 
reload(nmfb)
self = nmfb.NMFresult(source_dir)
self.set_output_dir(out=outdir)
self.load_NMF_H_matrix(Hfile)
self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
self.group_component_maps(match_field='pert_id')
# MI - 9 components
self.load_MI_matrix(MI_file_component)
self.MI_pairwise_comp(self.mi,match_field='pert_id')
self.intra_group_boxplot(space_name='c9')
self.boxplot_with_null(space_name='c9')
# MI - LM
self.load_MI_matrix(MI_file_inspace)
self.MI_pairwise_comp(self.mi,match_field='pert_id')
self.intra_group_boxplot(space_name='LM')
self.boxplot_with_null(space_name='LM')

################################################ 
### drug signatures projected in shRNA space ###
################################################ 

# output directory
outdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/VCAP_drug_c9_lm_epsilon/shRNA_projections/benchmark_graphs'
### Define Inputs
source_dir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/VCAP_drug_c9_lm_epsilon' # NMF result directory
nComponents = 9
prefix = 'clique_compound_classes'
dim = 'n380x978'
# Hfile = prefix +'_' + dim + '.H.k' + str(nComponents) + '.gct'
Hfile = prefix +'_' + dim + '.shRNA_projections.H.k' + str(nComponents) + '.gct'
WFile = prefix +'_' + dim + '.W.k' + str(nComponents) + '.gct'
MI_file_component = prefix + '_n448x978.MI.k' + str(nComponents) + '.gct'
MI_file_inspace = prefix + '_n448x978.MI.input_space.gct'
anntFile = prefix +'.v2.txt'
groupFile = '/xchip/cogs/projects/NMF/NMF_drug_shRNA/n69_drug_targets.gmt'

# run NMF module 
reload(nmfb)
self = nmfb.NMFresult(source_dir)
self.set_output_dir(out=outdir)
self.load_NMF_H_matrix(Hfile)
self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
self.group_component_maps(match_field='pert_id')
# MI - 9 components
self.load_MI_matrix(MI_file_component)
self.MI_pairwise_comp(self.mi,match_field='pert_id')
self.intra_group_boxplot(space_name='c9')
self.boxplot_with_null(space_name='c9')
# MI - LM
self.load_MI_matrix(MI_file_inspace)
self.MI_pairwise_comp(self.mi,match_field='pert_id')
self.intra_group_boxplot(space_name='LM')
self.boxplot_with_null(space_name='LM')

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
