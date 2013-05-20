#!/bin/sh

#run python dose anaysis tool
cellL=MCF7
timeP=24H
metric=wteslm

#plate control
inPut=/xchip/obelix/pod/brew/pc/ASG001_${cellL}_${timeP}/by_pert_id_pert_dose/ASG001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n85x978.gctx
outDir=/xchip/cogs/projects/ASG_dose_time/cmap_queries/$cellL/$timeP/$metric

#vehicle control
# inPut=/xchip/obelix/pod/brew/vc/ASG001_${cellL}_${timeP}/by_pert_id_pert_dose/ASG001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n85x978.gctx
# outDir=/xchip/cogs/projects/ASG_dose_time/cmap_queries/$cellL/$timeP/${metric}_vc

mkdir -p $outDir

#wteslm results
# rsltGctx=/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7/6H/wteslm/jan04/dose_analysis_tool.1357311873562/result_WTESLM.COMBINED_n85x398050.gctx
rsltGctx=/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7/24H/wteslm/jan04/dose_analysis_tool.1357324032334/result_WTESLM.COMBINED_n85x398050.gctx
# rsltGctx=/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/jan04/dose_analysis_tool.1357311906761/result_WTESLM.COMBINED_n85x398050.gctx
# rsltGctx=/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/24H/wteslm/jan04/dose_analysis_tool.1357311902816/result_WTESLM.COMBINED_n85x398050.gctx

# #run query + analysis
# python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_analysis_tool.py $inPut -o $outDir -q $metric -r
# # give existing query result gctx
python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_analysis_tool.py $inPut -o $outDir --resultDir $rsltGctx

#create a broad html page from the results
# ln -s /xchip/cogs/hogstrom/analysis/scratch/ASG_qqPlots/ example_lh
# ln -s /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/jan14/dose_analysis_tool.1358189363930/ PC3_6H_dose