'''
organize data for TA Pan cancer project

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
# wkdir = '/xchip/cogs/projects/NMF/TA_pan_cancer_OE_June_2014/TA_OE_ZSPCINF'
wkdir = '/xchip/cogs/projects/NMF/TA_pan_cancer_OE_May_2014/TA_OE_ZSPC_LM'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

processesed_type = 'ZSPC_LM' # 'COMPZ.MODZ_SCORE', , 'ZSPCINF'

#########################
### Run NMF projection ##
#########################

# COMPZ.MODZ_SCORE
nComponents = 20
dimDict = {'HA1E':'n3148x978'}

################
### load data ##
################

file_modz = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_COMPZ.MODZ_SCORE_n13974x22268.gctx'
file_qnorm = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_QNORM_n38534x978.gctx'
file_zspcinf = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_ZSPCINF_n38534x22268.gctx'

gt = gct.GCT(src=file_zspcinf)
gt.read()
ds = gt.frame

# signature annotations
sFile = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/inst.info'
sigInfo = pd.read_csv(sFile,sep='\t')
sigInfo.index = sigInfo.distil_id

####################
### prep sig ids ##
####################

### reduce to LM genes 
gLM = gct.GCT()
gLM.read_gctx_row_meta(src=file_qnorm) # load file with LM genes
lm_probes = gLM.get_rids()
ds_pan = ds.ix[ds.index.isin(lm_probes),:]

#split columns
colSer = pd.Series(ds_pan.columns)
colSer.name = 'sig_id'
colSplit = colSer.str.split('_')

# annotate acording to sig_id fields
colFrame = pd.DataFrame(colSer)
colFrame['plate'] = colSplit.apply(lambda x: x[0])
colFrame['cell_line'] = colSplit.apply(lambda x: x[1])
colFrame['tp'] = colSplit.apply(lambda x: x[2])
colFrame['rep'] = colSplit.apply(lambda x: x[3])
colFrame['well'] = colSplit.apply(lambda x: x[4])

plate_list = ['TA.OE009', 'TA.OE010', 'TA.OE011']
is_oe = colFrame.plate.isin(plate_list)
oe = colFrame[is_oe]

#####################
### WT update list ###
#####################

sheet_wt = '/xchip/cogs/projects/NMF/TA_pan_cancer_OE_May_2014/TA_OE_ZSPC_LM/HA1E/WT_update_analysis_Jul2014/PanCancerTemplateInfopDONR_partial.txt'
wtFrame = pd.read_csv(sheet_wt, sep='\t')
# filter out 'laters'
wtFrame = wtFrame[~(wtFrame['DerivedFrom(WTpDONR)'] == 'later')]
# annotate main gene 
splitSer = wtFrame.x_AllPanCancerMutations_YS.str.split("_")
mainGene = [x[0] for x in splitSer]
wtFrame['main_gene'] = mainGene
mainGrped = wtFrame.groupby('main_gene')
### get single DerivedFrom(WTpDONR) status for each gene
is_MUT = np.array(['>' in x for x in wtFrame.x_AllPanCancerMutations_YS])
mutFrm = wtFrame[is_MUT] 
formSet = set(mutFrm['DerivedFrom(WTpDONR)'])
# keep only WTs that match MUT open/closed form
is_match = wtFrame['DerivedFrom(WTpDONR)'].isin(list(formSet))
matchFrm = wtFrame[is_match]
matchFrm = matchFrm.sort('main_gene')

################################
### make gene signature gmt ###
################################
#x_allpancancermutations_ys
#x_annotgenesymbol
#x_template_gene_ys

# # reindex acording to OE plates
sigInfo = sigInfo.reindex(oe.sig_id)
sigInfo = sigInfo[sigInfo.pert_id.isin(matchFrm.pert_id)]
# sigGrped = sigInfo.groupby(['cell_id','pert_mfc_desc'])
cellGrped = sigInfo.groupby('cell_id')
for cellTup in cellGrped:
    cell = cellTup[0]
    cellFrm = cellTup[1]
    cellDir = wkdir + '/' + cell
    outF = cellDir + '/OE_annotations.txt'
    # reformat sig_id
    cellFrm['mod_sig_id'] = cellFrm.distil_id.str.replace(':','.')
    cellFrm.index = cellFrm.mod_sig_id
    cellFrm.to_csv(outF,sep='\t')
    ### make gene signature groups - gmt file
    # geneGrped = cellFrm.groupby('pert_mfc_desc')
    # geneGrped = cellFrm.groupby('x_mutation_status')
    # geneGrped = cellFrm.groupby('x_allpancancermutations_ys')
    geneGrped = cellFrm.groupby('pert_id')
    gmtList = []
    for grp in geneGrped:
        # leave out perturbations with less than 3 signatures
        if len(grp[1].index.values) > 2:
            gmtDictUp = {}
            gmtDictUp['desc'] = grp[0]
            # gmtDictUp['desc'] = str(list(set(grp[1].x_allpancancermutations_ys)))
            oc = grp[1]['x_openclosed_ys'].values[0]
            gn = list(set(grp[1].x_allpancancermutations_ys))[0]
            desc = gn + '_' + oc + '_' + grp[0]
            gmtDictUp['id'] = desc
            # gmtDictUp['desc'] = list(set(grp[1].x_allpancancermutations_ys))[0]
            gmtDictUp['sig'] = list(grp[1].index.values)
            gmtList.append(gmtDictUp)
    # gmtOut = cellDir + '/pert_id_oe_sig_id.gmt'
    gmtOut = cellDir + '/pert_id_open_closed_match.gmt'
    gmt.write(gmtList,gmtOut)

#########################
### run NMF benchmarks ##
#########################

# plate_dir_name = 'pert_id_grouped_analysis'
plate_dir_name = 'WT_update_analysis_Jul2014'
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
    groupFile = path1 + '/pert_id_open_closed_match.gmt'
    # groupFile = path1 + '/pert_id_oe_sig_id.gmt'
    # groupFile = path1 + '/mutation_status_oe_sig_id.gmt'
    # groupFile = path1 + '/' + grouping_file
    # run NMF module 
    reload(nmfb)
    self = nmfb.NMFresult(source_dir)
    self.set_output_dir(out=outdir)
    self.load_NMF_H_matrix(Hfile)
    self.load_NMF_W_matrix(WFile)
    self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
    self.load_input_matrix(prefix1+'.gct', modify_sig_id=True, reindex_by_H=True)
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
        has_WT = np.array(['>' not in x for x in grp_values])
        # select groups that contain WT and MUT
        if (sum(has_WT) > 0) & (sum(has_WT) < len(has_WT)):
            WT = grp_values[has_WT]
            MUT = grp_values[~has_WT]
            mutDict[gene] = tuple([WT,MUT])
    # # component heatmaps
    self.group_component_maps(match_field='signatures')
    ###################################
    ### pairwise comparisons of LM space
    ###################################
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
    self.intra_group_boxplot(space_name='LM_space',similarity_metric='mutual_information')
    self.boxplot_with_null(space_name='LM_space',similarity_metric='mutual_information')
    self.MUT_WT_comparison(self.pairwise_similarity_mtrx,mutDict,space_name='LM_space',
                        similarity_metric='rnkpt_MI',out_table=True)
    self.MUT_WT_boxplot(mutDict,space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='WT_MUT_boxplot_rnkpt_MI_LM_space',xlim_range=(-100,100),graph_title_str=prefix + ' - ')
    self.MUT_WT_graph_scatter(mutDict,space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_LM_space_wide',axis_scale=100,wt_median_thresh=None,
        graph_title_str=prefix + ' - ')
    self.MUT_WT_graph(mutDict,space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_LM_space_wide',axis_scale=100,wt_median_thresh=None,
        graph_title_str=prefix + ' - ')
    self.ONC_TSG_mapping(mutDict,space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_LM_space_wide',score_thresh=75)
    self.ONC_TSG_plot(space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='WT_MUT_diff_graphs_rnkpt_MI_LM_space_wide',axis_scale=100,
        graph_title_str=prefix + ' - ')
    self.ONC_TSG_ordered_plots(space_name='LM_space', similarity_metric='rnkpt_MI',
        out_graph_dir='connection_bins',xlim_range=(-100,100),axis_scale=100,
        graph_title_str=prefix + ' - ')
    ###################################
    ### perform CMAP queries with W space
    ###################################
    self.write_W_weights(n_probes_up=100,n_probes_dn=1,W_weight_dir='W_components')
    self.cmap_W_weight_query(W_weight_dir='W_components',up_file='W_component_high_weighted_probes.gmt',dn_file='W_component_low_weighted_probes.gmt')
    self.cmap_W_weight_summly(W_weight_dir='W_components')
