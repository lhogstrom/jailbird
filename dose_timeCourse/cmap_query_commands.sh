#!/bin/sh

# #query limited to is_gold
# #from: /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/24H_isgld
# rum -q local query_tool --uptag /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/24H_isgld/PC3_24H_signatures_LM_up.gmt --dntag /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/24H_isgld/PC3_24H_signatures_LM_dn.gmt --cid /xchip/cogs/projects/ASG_dose_time/cmap_queries/isGld_sig_ids.grp --metric eslm
# 
# #from: /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/6H_isgld
# rum -q local query_tool --uptag /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/6H_isgld/PC3_6H_signatures_LM_up.gmt --dntag /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/6H_isgld/PC3_6H_signatures_LM_dn.gmt --cid /xchip/cogs/projects/ASG_dose_time/cmap_queries/isGld_sig_ids.grp --metric eslm
# 
# 
# #from: /xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/6H_isgld
# rum -q local query_tool --uptag /xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/6H_isgld/MCF7_6H_signatures_LM_up.gmt --dntag /xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/6H_isgld/MCF7_6H_signatures_LM_dn.gmt --cid /xchip/cogs/projects/ASG_dose_time/cmap_queries/isGld_sig_ids.grp --metric eslm
# 
# #from: /xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/24H_isgld
# rum -q local query_tool --uptag /xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/24H_isgld/MCF7_24H_signatures_LM_up.gmt --dntag /xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/24H_isgld/MCF7_24H_signatures_LM_dn.gmt --cid /xchip/cogs/projects/ASG_dose_time/cmap_queries/isGld_sig_ids.grp --metric eslm


#run python dose anaysis tool
cellL=MCF7
timeP=6H
metric=wteslm

#plate control
inPut=/xchip/obelix/pod/brew/pc/ASG001_${cellL}_${timeP}/by_pert_id_pert_dose/ASG001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n85x978.gctx
outDir=/xchip/cogs/projects/ASG_dose_time/cmap_queries/$cellL/$timeP/$metric

#vehicle control
inPut=/xchip/obelix/pod/brew/vc/ASG001_${cellL}_${timeP}/by_pert_id_pert_dose/ASG001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n85x978.gctx
outDir=/xchip/cogs/projects/ASG_dose_time/cmap_queries/$cellL/$timeP/${metric}_vc

mkdir -p $outDir
#rsltGctx=/xchip/cogs/hogstrom/analysis/scratch/Nov20/dose_analysis_tool.1353449771597/nov20/my_analysis.query_tool.2012112017162991/result_ESLM.COMBINED_n85x398050.gctx
# rsltGctx=/xchip/cogs/projects/ASG_dose_time/cmap_queries/$cellL/$timeP/$metric/jan04/*/result_WTESLM.COMBINED_n85x398050.gctx

#wteslm results
# rsltGctx=/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7/6H/wteslm/jan04/dose_analysis_tool.1357311873562/result_WTESLM.COMBINED_n85x398050.gctx
# rsltGctx=/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7/24H/wteslm/jan04/dose_analysis_tool.1357324032334/result_WTESLM.COMBINED_n85x398050.gctx
# rsltGctx=/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/jan04/dose_analysis_tool.1357311906761/result_WTESLM.COMBINED_n85x398050.gctx
# rsltGctx=/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/24H/wteslm/jan04/dose_analysis_tool.1357311902816/result_WTESLM.COMBINED_n85x398050.gctx




# python /xchip/cogs/hogstrom/scripts/bptk3/bptk/cmap/tools/dose_analysis_tool.py $inPut -o $outDir -q $metric -rd $rsltGctx
# python /xchip/cogs/hogstrom/scripts/dose_timeCourse/ljh_dose_analysis_tool6.py  $inPut -o $outDir -q $metric -r
# python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_analysis_tool.py $inPut -o $outDir -q $metric -r
python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_analysis_tool.py $inPut -o $outDir --resultDir $rsltGctx

#create a broad html page from the results
# ln -s /xchip/cogs/hogstrom/analysis/scratch/ASG_qqPlots/ example_lh
# ln -s /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/jan14/dose_analysis_tool.1358189363930/ PC3_6H_dose