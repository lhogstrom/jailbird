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
wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_May_2014/TA_OE_ZSPCINF'
# wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_May_2014/TA_OE_ZSPC_LM'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

################
### load data ##
################

file_modz = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_COMPZ.MODZ_SCORE_n13968x22268.gctx'
file_qnorm = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_QNORM_n37799x978.gctx'
file_zspcinf = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_ZSPCINF_n37799x22268.gctx'

gt = gct.GCT(src=file_zspcinf)
gt.read()
ds = gt.frame

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



# signature subset 
file_lung_grp = '/cga/meyerson/brooks/TA/all_TA_for_jun10/all_TA_Lung_sig_ids.grp'
lungSigs = pd.read_csv(file_lung_grp)

# 
sFile = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/inst.info'
sigInfo = pd.read_csv(sFile,sep='\t')
sigInfo.index = sigInfo.distil_id

# important fileds in inst.info:
# pert_mfc_desc
# x_mutation_status
# x_preferredgenename
# x_tomconstructname

# # reindex acording to OE plates
# sigInfo = sigInfo.reindex(oe.sig_id)
# # sigGrped = sigInfo.groupby(['cell_id','pert_mfc_desc'])
# cellGrped = sigInfo.groupby('cell_id')
# for cellTup in cellGrped:
#     cell = cellTup[0]
#     cellFrm = cellTup[1]
#     cellDir = wkdir + '/' + cell
#     outF = cellDir + '/OE_annotations.txt'
#     # reformat sig_id
#     cellFrm['mod_sig_id'] = cellFrm.distil_id.str.replace(':','.')
#     cellFrm.index = cellFrm.mod_sig_id
#     cellFrm.to_csv(outF,sep='\t')
#     ### make gene signature groups - gmt file
#     # geneGrped = cellFrm.groupby('pert_mfc_desc')
#     geneGrped = cellFrm.groupby('x_mutation_status')
#     gmtList = []
#     for grp in geneGrped:
#         gmtDictUp = {}
#         gmtDictUp['id'] = grp[0]
#         # gmtDictUp['desc'] = grp[0]
#         gmtDictUp['desc'] = str(list(set(grp[1].x_mutation_status)))
#         gmtDictUp['sig'] = list(grp[1].index.values)
#         gmtList.append(gmtDictUp)
#     gmtOut = cellDir + '/mutation_status_oe_sig_id.gmt'
#     gmt.write(gmtList,gmtOut)


