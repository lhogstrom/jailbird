#!/bin/sh

#from matlab: 
# fname = '/xchip/obelix/pod/brew/pc/PRISM001_MCF7_24H/by_pert_id_pert_dose/PRISM001_MCF7_24H_COMPZ.MODZ_SCORE_LM_n59x978.gctx'
# db = parse_gctx(fname)
# rnk = rankorder(db.mat, 'direc', 'descend', 'fixties', false);
# rnkStruct = db;
# rnkStruct.mat = rnk;
# fout = '/xchip/cogs/hogstrom/analysis/scratch/prism/MCF7/24H/wteslm/feb05/dose_analysis_tool.1360084956001/PRISM001_MCF7_24H_rnk_LM.gctx'
# mkgctx(fout,rnkStruct);

rum -q local query_tool --score /xchip/obelix/pod/brew/pc/PRISM001_MCF7_24H/by_pert_id_pert_dose/PRISM001_MCF7_24H_COMPZ.MODZ_SCORE_LM_n59x978.gctx --rank_lm /xchip/cogs/hogstrom/analysis/scratch/prism/MCF7/24H/wteslm/feb05/dose_analysis_tool.1360084956001/PRISM001_MCF7_24H_rnk_LM_n59x978.gctx --uptag up_list.gmt --dntag dn_list.gmt --metric wteslm