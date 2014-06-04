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

################
### load data ##
################

file_modz = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_COMPZ.MODZ_SCORE_n13974x22268.gctx'
file_qnorm = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_QNORM_n38534x978.gctx'
file_zspcinf = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_ZSPCINF_n38534x22268.gctx'

gt = gct.GCT(src=file_zspcinf)
gt.read()
ds = gt.frame

# signature subset 
# file_lung_grp = '/cga/meyerson/brooks/TA/all_TA_for_jun10/all_TA_Lung_sig_ids.grp'
file_lung_grp = '/xchip/cga_home/brooks/TA/all_TA_for_jun10/all_TA_Lung_distil_ids.grp'
lungSigs = pd.read_csv(file_lung_grp,header=None, names=['sig_id'])
ds_lung = ds.reindex(columns=lungSigs.sig_id.values)

# signature annotations
sFile = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/inst.info'
sigInfo = pd.read_csv(sFile,sep='\t')
sigInfo.index = sigInfo.distil_id

#####################################
### prep GCT matrices by cell line ##
#####################################

processesed_type = 'ZSPC_LM' # 'COMPZ.MODZ_SCORE', , 'ZSPCINF'
### reduce to LM genes 
gLM = gct.GCT()
gLM.read_gctx_row_meta(src=file_qnorm) # load file with LM genes
lm_probes = gLM.get_rids()
ds_lung = ds_lung.ix[ds_lung.index.isin(lm_probes),:]

#split columns
colSer = pd.Series(ds_lung.columns)
colSer.name = 'sig_id'
colSplit = colSer.str.split('_')

# annotate acording to sig_id fields
colFrame = pd.DataFrame(colSer)
colFrame['plate'] = colSplit.apply(lambda x: x[0])
colFrame['cell_line'] = colSplit.apply(lambda x: x[1])
colFrame['tp'] = colSplit.apply(lambda x: x[2])
colFrame['rep'] = colSplit.apply(lambda x: x[3])
colFrame['well'] = colSplit.apply(lambda x: x[4])

################################
### make cell line gct files ###
################################

# save matrix for each cell line in OE experiments
is_oe = colFrame.plate.str.match('TA.OE0')
oe = colFrame[is_oe]
cell_grped = oe.groupby('cell_line')
for grpT in cell_grped:
    cell = grpT[0]
    cellDir = wkdir + '/' + cell
    if not os.path.exists(cellDir):
        os.mkdir(cellDir)
    grp = grpT[1]
    sigs = grp.sig_id
    cellFrm = ds_lung.ix[:,sigs.values]
    nGt = gct.GCT()
    nGt.build_from_DataFrame(cellFrm)
    outF = cellDir + '/' + cell + '_TA_JUN10_' + processesed_type
    nGt.write(outF,mode='gctx')

# convert gctx to gct
### run 'use Java-1.7' before running script ###
for cell in cell_grped.groups.keys():
    print(cell)
    cellDir = wkdir + '/' + cell
    outGCT = cellDir + '/' + cell
    globRes = glob.glob(outGCT+'*.gctx')
    print(globRes[0])
    cmd2 = 'convert-dataset -i ' + globRes[0]
    os.system(cmd2)

################################
### make gene signature gmt ###
################################

# # reindex acording to OE plates
sigInfo = sigInfo.reindex(oe.sig_id)
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
    geneGrped = cellFrm.groupby('x_mutation_status')
    gmtList = []
    for grp in geneGrped:
        gmtDictUp = {}
        gmtDictUp['id'] = grp[0]
        # gmtDictUp['desc'] = grp[0]
        gmtDictUp['desc'] = str(list(set(grp[1].x_mutation_status)))
        gmtDictUp['sig'] = list(grp[1].index.values)
        gmtList.append(gmtDictUp)
    gmtOut = cellDir + '/mutation_status_oe_sig_id.gmt'
    gmt.write(gmtList,gmtOut)

#########################
### Run NMF projection ##
#########################

# COMPZ.MODZ_SCORE
nComponents = 20
dimDict = {'A549':'n4487x978', # 
'AALE':'n2235x978',
'H1299':'n1503x978',
'SALE':'n2128x978'}

# ZSPCINF
# nComponents = 20
# dimDict = {'A375':'n2245x22268',
# 'A549':'n5608x22268', # 
# 'AALE':'n3356x22268',
# 'H1299':'n2597x22268',
# 'HA1E':'n5371x22268',
# 'PC3':'n2246x22268',
# 'SALE':'n3215x22268'}

#specifications for subprocess
processes = set()
max_processes = 9 
### run jobs
for cell in cell_grped.groups.keys():
    print cell
    dim = dimDict[cell]
    arg1 = wkdir + '/' + cell # working directory
    # arg2 = cell + '_TA_JUN10_COMPZ.MODZ_SCORE_' + dim
    # arg2 = cell + '_TA_JUN10_ZSPCINF_' + dim
    arg2 = cell + '_TA_JUN10_' + processesed_type + '_' + dim
    cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/NMF_code_no_viz.v2.R', # 
         arg1,
         arg2])
    # os.system(cmd)
    processes.add(subprocess.Popen(cmd,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes if p.poll() is not None)

########################
### load annotations ###
########################

# load lung-driver sheet
sheet_driver = '/xchip/cogs/projects/NMF/TA_lung_OE_May_2014/LUAD_ORFs_UpdatedAnnotations_GoogleSheet-TA.OE012_013_ORF_GENE_MUT_CAT.tsv'
drivers = pd.read_csv(sheet_driver,sep='\t')
#x_mutation_category_20140122
cat_grped = drivers.groupby('x_mutation_category_20140122')
count_dict = {}
for grp in cat_grped:
    ptype = grp[0]
    mtx = grp[1]
    sig_match = sigInfo[sigInfo.pert_id.isin(mtx.pert_id)]
    cell_grp = sig_match.groupby('cell_id')
    if len(cell_grp.indices) > 0:
        count_dict[ptype] = cell_grp.apply(len)
cell_counts = pd.DataFrame(count_dict)
outF = os.path.join(wkdir,'mutation_category_counts.txt')
cell_counts.to_csv(outF,sep='\t')

# important fileds in inst.info:
# pert_mfc_desc
# x_mutation_status
# x_preferredgenename
# x_tomconstructname

#########################
### run NMF benchmarks ##
#########################

for prefix in dimDict:
    dim = dimDict[prefix]
    path1 = wkdir + '/' + prefix
    prefix1 = prefix + '_TA_JUN10_'+ processesed_type + '_' + dim
    outdir = path1 + '/mutation_status_benchmark_graphs'
    source_dir = path1
    Hfile = prefix1 + '.H.k' + str(nComponents) + '.gct'
    # WFile = prefix1 + '.W.k' + str(nComponents) + '.gct'
    MI_file_component = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.k' + str(nComponents) + '.gct'
    MI_file_inspace = prefix + '_TA_JUN10_'+ processesed_type + '_n.MI.input_space.gct'
    anntFile = 'OE_annotations.txt'
    # groupFile = path1 + '/gene_oe_sig_id.gmt'
    groupFile = path1 + '/mutation_status_oe_sig_id.gmt'
    # run NMF module 
    reload(nmfb)
    self = nmfb.NMFresult(source_dir)
    self.set_output_dir(out=outdir)
    self.load_NMF_H_matrix(Hfile)
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
        has_WT = np.array(['WT' in x for x in grp_values])
        # select groups that contain WT and MUT
        if (sum(has_WT) > 0) & (sum(has_WT) < len(has_WT)):
            WT = grp_values[has_WT]
            MUT = grp_values[~has_WT]
            mutDict[gene] = tuple([WT,MUT])
    ###################################
    ### pairwise comparisons of H space
    ###################################
    # Pearson Corr 
    self.calculate_corr_matrix(H_mtrx_corr=True) # pairwise corr on H-matrix
    self.MI_pairwise_comp(self.pairwise_corr,match_field='signatures',out_table=True)
    # self.intra_group_boxplot(space_name='20_components',similarity_metric='Pearson correlation')
    # self.boxplot_with_null(space_name='20_components',similarity_metric='Pearson correlation')
    self.MUT_WT_comparison(self.pairwise_corr,mutDict,space_name='20_components',
                        similarity_metric='Pearson_correlation',out_table=True)
    self.MUT_WT_graph(mutDict,space_name='input_space', similarity_metric='Pearson_correlation',out_graph_dir='WT_MUT_graphs_Pearson_c20')
    # Mutual information
    self.load_MI_matrix(MI_file_component)
    self.MI_pairwise_comp(self.mi,match_field='signatures',out_table=True)
    self.intra_group_boxplot(space_name='20_components',similarity_metric='mutual_information')
    self.boxplot_with_null(space_name='20_components',similarity_metric='mutual_information')
    self.MUT_WT_comparison(self.mi,mutDict,space_name='20_components',
                        similarity_metric='mutual_information',out_table=True)
    self.MUT_WT_graph(mutDict,space_name='20_components', similarity_metric='mutual_information',out_graph_dir='WT_MUT_graphs_MI_c20')
    ###################################
    ### pairwise comparisons of LM space
    ###################################
    self.calculate_corr_matrix(H_mtrx_corr=False) # pairwise corr on input matrix
    self.MI_pairwise_comp(self.pairwise_corr,match_field='signatures')
    self.intra_group_boxplot(space_name='LM_space',similarity_metric='Pearson_correlation')
    self.boxplot_with_null(space_name='LM_space',similarity_metric='Pearson_correlation')
    # component heatmaps
    # self.group_component_maps(match_field='signatures')
    self.MUT_WT_comparison(self.pairwise_corr,mutDict,space_name='LM_space',
                        similarity_metric='Pearson_correlation',out_table=True)
    self.MUT_WT_graph(mutDict,space_name='LM_space', similarity_metric='Pearson_correlation',out_graph_dir='WT_MUT_graphs_Pearson_LM')
    # Mutual information
    self.load_MI_matrix(MI_file_inspace)
    self.MI_pairwise_comp(self.mi,match_field='signatures',out_table=True)
    self.intra_group_boxplot(space_name='LM_space',similarity_metric='mutual_information')
    self.boxplot_with_null(space_name='LM_space',similarity_metric='mutual_information')
    self.MUT_WT_comparison(self.mi,mutDict,space_name='LM_space',
                        similarity_metric='mutual_information',out_table=True)
    self.MUT_WT_graph(mutDict,space_name='LM_space', similarity_metric='mutual_information',
        out_graph_dir='WT_MUT_graphs_MI_LM_space',WT_MUT_comparison=True,wt_median_thresh=.3)
    self.MUT_WT_graph(mutDict,space_name='LM_space', similarity_metric='mutual_information',
        out_graph_dir='WT_MUT_graphs_MI_LM_space',WT_MUT_comparison=False,wt_median_thresh=.3)

