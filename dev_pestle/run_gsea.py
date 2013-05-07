#!/bin/py
'''
run gsea using the python tool
'''
import cmap.io.gct as gct
import cmap.analytics.gsea as gsea
import matplotlib.pyplot as plt
import numpy as np
import os

data1 = '/xchip/cogs/projects/PRISM/dose_plate_output-by_pert_id_pert_dose/MCF7/24H/vc/apr17/dose_plate_tool.1366234068487/CMAP-AZD-1152HQPA_linear_heatmap_data.gctx'
#LM
gctLM = '/xchip/cogs/data/brew/a2y13q1/PRISM001_PC3_6H/by_pert_id_pert_dose/PRISM001_PC3_6H_COMPZ.MODZ_SCORE_LM_n59x978.gctx'
#INF
gctINF = '/xchip/cogs/data/brew/a2y13q1/PRISM001_PC3_6H/by_pert_id_pert_dose/PRISM001_PC3_6H_COMPZ.MODZ_SCORE_n59x22268.gctx'
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/gsea_tests'

#need rank ordered list of probes with z-scores, then pre-ranked
db = gct.GCT()
# db.read(gctINF)
db.read_gctx_col_meta(gctINF)
db.read_gctx_matrix(gctINF)
probes = db.get_gctx_rid()
# db.read_gctx_row_meta(gctINF)
#order probes acording to a single perturbation - store z-score
# 'GSK-1070916'
db.get_column_meta('pert_dose')[34]
zCol = db.matrix[:,34]
iSortZ = np.argsort(zCol)[::-1]
zSort = zCol[iSortZ]
sProbes = [probes[i] for i in iSortZ]
#write probes to file
fout = os.path.join(work_dir,'probes_gsk_1070916_10um.txt')
with open(fout,'w') as f:
	for i,probe in enumerate(sProbes):
		z = zSort[i]
		f.write(probe + '\t')
		f.write(str(z) + '\t')
		f.write('\n')


