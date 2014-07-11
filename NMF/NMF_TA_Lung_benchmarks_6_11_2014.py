'''
organize data for TA lung cancer project

'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
import subprocess
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import cmap.io.gmt as gmt
import cmap.analytics.NMF_benchmarks as nmfb

######################
### set working dir ##
######################

# wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_May_2014/TA_OE_qnorm'
wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_June_2014/TA_OE_ZSPCINF'
# wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_May_2014/TA_OE_ZSPC_LM'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

processesed_type = 'ZSPC_LM' # 'COMPZ.MODZ_SCORE', , 'ZSPCINF'

#########################
### Run NMF projection ##
#########################

# COMPZ.MODZ_SCORE
nComponents = 20
dimDict = {'A549':'n4487x978', # 
'AALE':'n2235x978',
'H1299':'n1503x978',
'SALE':'n2128x978'}

#############################
### MI to rnkpt conversion ##
#############################

# for prefix in dimDict:
#     print prefix
#     dim = dimDict[prefix]
#     path1 = wkdir + '/' + prefix
#     prefix1 = prefix + '_TA_JUN10_'+ processesed_type + '_' + dim
#     source_dir = path1
#     Hfile = prefix1 + '.H.k' + str(nComponents) + '.gct'
#     # WFile = prefix1 + '.W.k' + str(nComponents) + '.gct'
#     MI_file_component = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.k' + str(nComponents) + '.gct'
#     MI_file_inspace = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.input_space.gct'
#     ### convert gct v.2 to gctx
#     reload(nmfb)
#     self = nmfb.NMFresult(source_dir)
#     self.load_MI_matrix(MI_file_component,gctx_out_file=prefix+'_TA_JUN10_ZSPC_LM_n.MI.k20')
#     self.load_MI_matrix(MI_file_inspace,gctx_out_file=prefix+'_TA_JUN10_ZSPC_LM_n.MI.input_space')
#     ### convert MI matrix to MI_rnkpt
#     MI_gctx_inspace = glob.glob(path1 + '/*MI.input_space*.gctx')[0]
#     MI_gctx_component = glob.glob(path1 + '/*MI.k*.gctx')[0]
#     rnkpt_file_component = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.rnkpt.k' + str(nComponents) + '.gctx'
#     rnkpt_file_inspace = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.rnkpt.input_space.gctx'
#     cmd_pre = 'python /xchip/cogs/hogstrom/scripts/jailbird/NMF/transform2rankpoint.py -i '
#     cmd_str_1 = cmd_pre + os.path.join(path1, MI_gctx_inspace) + ' -o ' + os.path.join(path1, rnkpt_file_inspace)
#     cmd_str_2 = cmd_pre + os.path.join(path1, MI_gctx_component) + ' -o ' + os.path.join(path1, rnkpt_file_component)
#     os.system(cmd_str_1)
#     os.system(cmd_str_2)

#########################
### reindex by plate ##
#########################

## 'TA.OE012','TA.OE013' plates
plate_match = ['TA.OE012','TA.OE013']
plate_dir_name = 'Plate_12_13_analysis'
grouping_file = 'mutation_status_plate_12_13_oe_sig_id.gmt'
# 'TA.OE006' plates
# plate_match = ['TA.OE006']
# plate_dir_name = 'Plate_006_analysis'
# grouping_file = 'mutation_status_plate_006_oe_sig_id.gmt'
# dimDict = {'A549':'n4487x978'}

# signatures without enough replicates
old_lung_grp = '/xchip/cga_home/brooks/TA/all_TA_for_jun10/depracated/all_TA_Lung_distil_ids.grp'
oldSigs = pd.read_csv(old_lung_grp,header=None, names=['sig_id'])
new_lung_grp = '/cga/meyerson/brooks/TA/all_TA_for_jun10/all_TA_Lung_distil_ids.grp'
newSigs = pd.read_csv(new_lung_grp,header=None, names=['sig_id'])
# get only the new sigs which appear in the old list
oldSet = set(oldSigs.sig_id.values)
newSet = set(newSigs.sig_id.values)
fullSet = oldSet.intersection(newSet)
lungSigs = pd.Series(list(fullSet))
lungSigs.name = 'sig_id'
colSplit = lungSigs.str.split('_')

# annotate acording to sig_id fields
colFrame = pd.DataFrame(lungSigs)
colFrame['plate'] = colSplit.apply(lambda x: x[0])
colFrame['cell_line'] = colSplit.apply(lambda x: x[1])
colFrame['tp'] = colSplit.apply(lambda x: x[2])
colFrame['rep'] = colSplit.apply(lambda x: x[3])
colFrame['well'] = colSplit.apply(lambda x: x[4])
# select plates
matchFrm = colFrame[colFrame.plate.isin(plate_match)]
matchGrp = matchFrm.groupby('cell_line')

################################
### make gene signature gmt ###
################################
# signature annotations
sFile = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/inst.info'
sigInfo = pd.read_csv(sFile,sep='\t')
sigInfo.index = sigInfo.distil_id

# # reindex acording to OE plates
sigInfo = sigInfo.reindex(matchFrm.sig_id)
# sigGrped = sigInfo.groupby(['cell_id','pert_mfc_desc'])
cellGrped = sigInfo.groupby('cell_id')
for cellTup in cellGrped:
    cell = cellTup[0]
    cellFrm = cellTup[1]
    cellFrm['mod_sig_id'] = cellFrm.distil_id.str.replace(':','.')
    cellFrm.index = cellFrm.mod_sig_id
    cellDir = wkdir + '/' + cell
    ### make gene signature groups - gmt file
    geneGrped = cellFrm.groupby('x_mutation_status')
    gmtList = []
    for grp in geneGrped:
        gmtDictUp = {}
        gmtDictUp['id'] = grp[0]
        # gmtDictUp['desc'] = grp[0]
        gmtDictUp['desc'] = str(list(set(grp[1].x_mutation_status)))
        gmtDictUp['sig'] = list(grp[1].index.values)
        gmtList.append(gmtDictUp)
    gmtOut = cellDir + '/' + grouping_file
    gmt.write(gmtList,gmtOut)

#########################
### run NMF benchmarks ##
#########################

for prefix in dimDict:
    print prefix
    dim = dimDict[prefix]
    path1 = wkdir + '/' + prefix
    prefix1 = prefix + '_TA_JUN10_'+ processesed_type + '_' + dim
    # outdir = path1 + '/mutation_status_benchmark_graphs'
    outdir = path1 + '/' + plate_dir_name
    source_dir = path1
    Hfile = prefix1 + '.H.k' + str(nComponents) + '.gct'
    WFile = prefix1 + '.W.k' + str(nComponents) + '.gct'
    MI_file_component = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.k' + str(nComponents) + '.gct'
    MI_file_inspace = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.input_space.gct'
    MI_rnkpt_component = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.rnkpt.k' + str(nComponents) + '.gctx'
    MI_rnkpt_inspace = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.rnkpt.input_space.gctx'
    anntFile = 'OE_annotations.txt'
    # groupFile = path1 + '/gene_oe_sig_id.gmt'
    # groupFile = path1 + '/mutation_status_oe_sig_id.gmt'
    groupFile = path1 + '/' + grouping_file
    # run NMF module 
    reload(nmfb)
    self = nmfb.NMFresult(source_dir)
    self.set_output_dir(out=outdir)
    self.load_NMF_H_matrix(Hfile)
    self.load_NMF_W_matrix(WFile)
    self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
    self.load_input_matrix(prefix1+'.gct', modify_sig_id=True, reindex_by_H=True)
    # re-write W matrix with gene symbols 
    W_symbol_out = prefix + '_TA_JUN10_'+ processesed_type + '.gene_symbols.W.k' + str(nComponents)
    self.probe_id_to_gene_symbol(self.Wmtrx,outfile=W_symbol_out)
    #########################
    ### WT vs MUT groupings ##
    #########################
    colSplit = self.groupFrm.id.str.split("_")
    # annotate acording to sig_id fields
    colFrame = pd.DataFrame(self.groupFrm.id)
    colFrame['main_gene'] = colSplit.apply(lambda x: x[0])
    colFrame.index = self.groupFrm.id
    geneGrped = colFrame.groupby('main_gene')
    mutDict = {}
    for grp_tup in geneGrped:
        gene = grp_tup[0]
        grp = grp_tup[1]
        grp_values = grp.id.values
        is_WT2 = np.array(['WT.2' in x for x in grp_values])
        grp_values = grp_values[~is_WT2] # remove WT2
        # has_WT = np.array(['WT' in x for x in grp_values])
        has_WT = np.array(['WT' in x for x in grp_values])
        # select groups that contain WT and MUT
        if (sum(has_WT) > 0) & (sum(has_WT) < len(has_WT)):
            WT = grp_values[has_WT]
            MUT = grp_values[~has_WT]
            mutDict[gene] = tuple([WT,MUT])
    # # component heatmaps
    # self.group_component_maps(match_field='signatures')
    ###################################
    ### pairwise comparisons of H space
    ###################################
    # Pearson Corr 
    # self.calculate_corr_matrix(H_mtrx_corr=True) # pairwise corr on H-matrix
    # self.MI_pairwise_comp(self.pairwise_corr,match_field='signatures',out_table=True)
    # # self.intra_group_boxplot(space_name='20_components',similarity_metric='Pearson correlation')
    # # self.boxplot_with_null(space_name='20_components',similarity_metric='Pearson correlation')
    # self.MUT_WT_comparison(self.pairwise_corr,mutDict,space_name='20_components',
    #                     similarity_metric='Pearson_correlation',out_table=True)
    # self.MUT_WT_graph(mutDict,space_name='input_space', similarity_metric='Pearson_correlation',out_graph_dir='WT_MUT_graphs_Pearson_c20')
    # Mutual information
    # self.load_MI_matrix(MI_file_component)
    # self.MI_pairwise_comp(self.mi,match_field='signatures',out_table=True)
    # self.intra_group_boxplot(space_name='20_components',similarity_metric='mutual_information')
    # self.boxplot_with_null(space_name='20_components',similarity_metric='mutual_information')
    # self.MUT_WT_comparison(self.mi,mutDict,space_name='20_components',
    #                     similarity_metric='mutual_information',out_table=True)
    # self.MUT_WT_graph_scatter(mutDict,space_name='20_components', similarity_metric='mutual_information',
    #     out_graph_dir='WT_MUT_diff_graphs_MI_c20',wt_median_thresh=None,graph_title_str=prefix + ' - ')
    # self.MUT_WT_graph(mutDict,space_name='20_components', similarity_metric='mutual_information',
    #     out_graph_dir='WT_MUT_graphs_MI_c20',WT_MUT_comparison=False,wt_median_thresh=.4)
    # Mutual information - rankpoint 
    # self.load_similarity_matrix(MI_rnkpt_component, similarity_metric='rnkpt_MI')
    # self.MI_pairwise_comp(self.pairwise_similarity_mtrx,match_field='signatures',out_table=True)
    # self.intra_group_boxplot(space_name='20_components',similarity_metric='rnkpt_MI',
    #         xlim_range=(-100,100))
    # self.boxplot_with_null(space_name='20_components',similarity_metric='rnkpt_MI')
    # self.MUT_WT_comparison(self.pairwise_similarity_mtrx,mutDict,space_name='20_components',
    #                     similarity_metric='rnkpt_MI',out_table=True)
    # self.MUT_WT_boxplot(mutDict,space_name='20_components', similarity_metric='rnkpt_MI',
    #     out_graph_dir='WT_MUT_boxplot_rnkpt_MI_c20',xlim_range=(-100,100),graph_title_str=prefix + ' - ')
    # self.MUT_WT_graph_scatter(mutDict,space_name='20_components', similarity_metric='rnkpt_MI',
    #     out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_c20',axis_scale=100,wt_median_thresh=None,
    #     graph_title_str=prefix + ' - ')
    # self.MUT_WT_graph(mutDict,space_name='20_components', similarity_metric='rnkpt_MI',
    #     out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_c20',axis_scale=100,wt_median_thresh=None,
    #     graph_title_str=prefix + ' - ')
    ###################################
    ### pairwise comparisons of LM space
    ###################################
    # Pearson Corr 
    # self.calculate_corr_matrix(H_mtrx_corr=False) # pairwise corr on input matrix
    # self.MI_pairwise_comp(self.pairwise_corr,match_field='signatures')
    # self.intra_group_boxplot(space_name='LM_space',similarity_metric='Pearson_correlation')
    # self.boxplot_with_null(space_name='LM_space',similarity_metric='Pearson_correlation')
    # self.MUT_WT_comparison(self.pairwise_corr,mutDict,space_name='LM_space',
    #                     similarity_metric='Pearson_correlation',out_table=True)
    # self.MUT_WT_graph(mutDict,space_name='LM_space', similarity_metric='Pearson_correlation',out_graph_dir='WT_MUT_graphs_Pearson_LM')
    # Mutual information
    # self.load_MI_matrix(MI_file_inspace)
    # self.MI_pairwise_comp(self.mi,match_field='signatures',out_table=True)
    # self.intra_group_boxplot(space_name='LM_space',similarity_metric='mutual_information')
    # self.boxplot_with_null(space_name='LM_space',similarity_metric='mutual_information')
    # self.MUT_WT_comparison(self.mi,mutDict,space_name='LM_space',
    #                     similarity_metric='mutual_information',out_table=True)
    # self.MUT_WT_graph(mutDict,space_name='LM_space', similarity_metric='mutual_information',
    #     out_graph_dir='WT_MUT_diff_graphs_MI_LM_space',wt_median_thresh=None,graph_title_str=prefix + ' - ')
    # self.MUT_WT_graph_scatter(mutDict,space_name='LM_space', similarity_metric='mutual_information',
    #     out_graph_dir='WT_MUT_diff_graphs_MI_LM_space',wt_median_thresh=None,graph_title_str=prefix + ' - ')
    # self.MUT_WT_boxplot(mutDict,space_name='LM_space', similarity_metric='mutual_information',
    #     out_graph_dir='WT_MUT_boxplot_MI_LM_space',graph_title_str=prefix + ' - ')
    # Mutual information - rankpoint 
    self.load_similarity_matrix(MI_rnkpt_inspace, similarity_metric='rnkpt_MI',reindex_ids=None)
    self.MI_pairwise_comp(self.pairwise_similarity_mtrx,match_field='signatures',out_table=True)
    # self.intra_group_boxplot(space_name='LM_space',similarity_metric='mutual_information')
    # self.boxplot_with_null(space_name='LM_space',similarity_metric='mutual_information')
    self.MUT_WT_comparison(self.pairwise_similarity_mtrx,mutDict,space_name='LM_space',
                        similarity_metric='rnkpt_MI',out_table=True)
    # self.MUT_WT_boxplot(mutDict,space_name='LM_space', similarity_metric='rnkpt_MI',
    #     out_graph_dir='WT_MUT_boxplot_rnkpt_MI_LM_space',xlim_range=(-100,100),graph_title_str=prefix + ' - ')
    self.MUT_WT_graph_scatter(mutDict,space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_LM_space',axis_scale=100,wt_median_thresh=None,
        graph_title_str=prefix + ' - ')
    self.MUT_WT_graph(mutDict,space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_LM_space_wide',axis_scale=100,wt_median_thresh=None,
        graph_title_str=prefix + ' - ')
    self.ONC_TSG_mapping(mutDict,space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_LM_space_wide',score_thresh=90)
    self.ONC_TSG_plot(space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_LM_space_wide',axis_scale=100,
        graph_title_str=prefix + ' - ')
    # self.ONC_TSG_ordered_plots(space_name='LM_space', similarity_metric='rnkpt_MI',
    #     out_graph_dir='connection_bins',xlim_range=(-100,100),axis_scale=100,
    #     graph_title_str=prefix + ' - ')
    ###################################
    ### perform CMAP queries with W space
    ###################################
    self.write_W_weights(n_probes_up=100,n_probes_dn=1,W_weight_dir='W_components')
    self.cmap_W_weight_query(W_weight_dir='W_components',up_file='W_component_high_weighted_probes.gmt',dn_file='W_component_low_weighted_probes.gmt')
    self.cmap_W_weight_summly(W_weight_dir='W_components')

### make hyperlinks for query results
for prefix in dimDict:
    print prefix
    dim = dimDict[prefix]
    path1 = wkdir + '/' + prefix
    prefix1 = prefix + '_TA_JUN10_'+ processesed_type + '_' + dim
    # outdir = path1 + '/mutation_status_benchmark_graphs'
    outdir = path1 + '/' + plate_dir_name + '/W_components'
    outpath = glob.glob('/'.join([outdir,'jul02/*']))[0]
    hyperLnkPath = '/xchip/cogs/web/icmap/hogstrom/TA_Lung_component_annotation/Plate_12_13_analysis/' + prefix
    cmd = ' '.join(['ln -s',
             outpath,
             hyperLnkPath])
    os.system(cmd)

