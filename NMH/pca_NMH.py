#!/bin/py


import cmap.analytics.pca as pca
import cmap.io.gct as gct

NMH001_NEU_6H = gct.GCT('/xchip/obelix/pod/brew/vc/NMH001_NEU_6H/by_pert_id/NMH001_NEU_6H_COMPZ.MODZ_SCORE_LM_n336x978.gctx',read=True)
pca.gct_pca(NMH001_NEU_6H)



#cluster correlation matrix from GENE-E
fnameClust = '/xchip/cogs/projects/NMH/cfwork/NMH001_NEU_6H_spearman_corr_matrix_clustered.gct'
clust_NMH001_NEU_6H = gct.GCT(fnameClust)