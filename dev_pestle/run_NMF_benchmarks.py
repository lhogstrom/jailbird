#! /usr/bin/env python

'''
Example script to run NMF_benchmarks module

LH 02/2014
'''
import os
import cmap.analytics.NMF_benchmarks as nmfb

###############################
### example 1 - Drug Groups ###
###############################

# output directory
outdir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation2/NMF_benchmark_development_module/A375_c9_lm_epsilon'
### Define Inputs
source_dir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation2/A375_c9_lm_epsilon' # NMF result directory
nComponents = 9
prefix = 'A375_c9_lm_epsilon'
dim = 'n473x978'
Hfile = 'clique_compound_classes_' + dim + '.H.k' + str(nComponents) + '.gct'
WFile = 'clique_compound_classes_' + dim + '.W.k' + str(nComponents) + '.gct'
MI_file_component = 'clique_compound_classes.MI.k' + str(nComponents) + '.gct'
MI_file_inspace = 'clique_compound_classes.MI.input_space.gct'
anntFile = 'clique_compound_classes.v2.txt'
groupFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140221/cliques.gmt'

# run NMF module 
reload(nmfb)
self = nmfb.NMFresult(source_dir)
self.set_output_dir(out=outdir)
self.load_NMF_H_matrix(Hfile)
self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
# self.group_component_maps(match_field='pert_id')
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

##############################
### example 2 - TA lung OE ###
##############################

# output directory
outdir = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/oe_pilot/NMF_projection_c9/module_benchmarks'
### Define Inputs
source_dir = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/oe_pilot/NMF_projection_c9'
# prefix = 'NMF_projection_c9'
nComponents = 9
modzPrefix = 'TA.OE006_A549_96H_mod_sig_id'
Hfile = modzPrefix + '.H.k' + str(nComponents) + '.gct'
WFile = modzPrefix + '.W.k' + str(nComponents) + '.gct'
MI_file_component = modzPrefix + '.MI.k' + str(nComponents) + '.gct'
MI_file_inspace = modzPrefix + '.MI.input_space.gct'
anntFile = 'TA.OE006_A549_96H_mod_sig_id_annotations2.txt'
groupFile = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/oe_pilot/core_lung_drivers_sig_id.gmt'

# run NMF module 
reload(nmfb)
self = nmfb.NMFresult(source_dir)
self.set_output_dir(out=outdir)
self.load_NMF_H_matrix(Hfile)
self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
# self.group_component_maps(match_field='signatures')
# MI - 9 components
self.load_MI_matrix(MI_file_component)
self.MI_pairwise_comp(self.mi,match_field='signatures')
self.intra_group_boxplot(space_name='c9')
self.boxplot_with_null(space_name='c9')
# MI - LM
self.load_MI_matrix(MI_file_inspace)
self.MI_pairwise_comp(self.mi,match_field='signatures')
self.intra_group_boxplot(space_name='LM')
self.boxplot_with_null(space_name='LM')

#################################
### example 3 - TA lung shRNA ###
#################################

# output directory
outdir = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/shRNA_pilot/NMF_projection_c20/module_benchmarks'
### Define Inputs
source_dir = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/shRNA_pilot/NMF_projection_c20'
nComponents = 20
modzPrefix = 'TA.KD009_KD010_A549_96H_mod_sig_id'
Hfile = modzPrefix + '.H.k' + str(nComponents) + '.gct'
WFile = modzPrefix + '.W.k' + str(nComponents) + '.gct'
MI_file_component = modzPrefix + '.MI.k' + str(nComponents) + '.gct'
MI_file_inspace = modzPrefix + '.MI.input_space.gct'
anntFile = 'TA.KD009_KD010_A549_96H_annotations.txt'
groupFile = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/shRNA_pilot/gene_shRNA_sig_id.gmt'

# run NMF module 
reload(nmfb)
self = nmfb.NMFresult(source_dir)
self.set_output_dir(out=outdir)
self.load_NMF_H_matrix(Hfile)
self.load_annotations(anntFile,sig_col=1,drop_extra_signatures=True,signature_group_file=groupFile)
# self.group_component_maps(match_field='signatures')
# MI - components
self.load_MI_matrix(MI_file_component)
self.MI_pairwise_comp(self.mi,match_field='signatures')
self.intra_group_boxplot(space_name='c9')
self.boxplot_with_null(space_name='c9')
# MI - LM
self.load_MI_matrix(MI_file_inspace)
self.MI_pairwise_comp(self.mi,match_field='signatures')
self.intra_group_boxplot(space_name='LM')
self.boxplot_with_null(space_name='LM')

#########################################
### example 4 - shRNA of drug targets ###
#########################################

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


