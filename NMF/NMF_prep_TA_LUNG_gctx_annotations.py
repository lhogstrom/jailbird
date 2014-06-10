'''
re-write 

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

### cell line gcts w/ annotation
# gFile = '/xchip/cogs/web/icmap/custom/TA/brew/pc/TA.OE013_A549_96H/TA.OE013_A549_96H_QNORM_n1117x978.gctx'
# gt_plate = gct.GCT(src=gFile)
# gt_plate.read()
# ds_plate = gt_plate.frame

file_modz = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_COMPZ.MODZ_SCORE_n13974x22268.gctx'
file_qnorm = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_QNORM_n38534x978.gctx'
file_zspcinf = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_ZSPCINF_n38534x22268.gctx'

gt = gct.GCT(src=file_modz)
gt.read()
ds = gt.frame

wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_June_2014/gctx_files_annotated/MODZ_INF'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

# # save matrix for each cell line in OE experiments
# is_oe = colFrame.plate.str.match('TA.OE0')
# oe = colFrame[is_oe]
# cell_grped = oe.groupby('cell_line')
# for grpT in cell_grped:
#     cell = grpT[0]
#     cellDir = wkdir + '/' + cell
#     if not os.path.exists(cellDir):
#         os.mkdir(cellDir)
#     grp = grpT[1]
#     sigs = grp.sig_id
#     nGt = gct.GCT(src=file_modz)
#     nGt.read(cid=sigs)
#     # ds = gt.frame
#     outF = cellDir + '/' + cell + '_TA_JUN10_' + processesed_type
#     nGt.write(outF,mode='gctx')

processed_type = 'MODZ_INF'
cell_lines = ['A549', 'H1299', 'AALE', 'SALE']
for cell in cell_lines:
    cellDir = wkdir + '/' + cell
    if not os.path.exists(cellDir):
        os.mkdir(cellDir)
    sFile = '/cga/meyerson/brooks/TA/all_TA_for_jun10/sig_ids_by_gene_w_OE007/all_TA_Lung_w_OE007_' + cell + '_sig_ids.grp'
    sigs = pd.read_csv(sFile,header=None)
    nGt = gct.GCT(src=file_modz)
    # nGt.read(cid=list(sigs[0].values),rid='lm_epsilon')
    nGt.read(cid=list(sigs[0].values))
    outF = cellDir + '/' + cell + '_TA_JUN10_' + processesed_type
    nGt.write(outF,mode='gctx')
