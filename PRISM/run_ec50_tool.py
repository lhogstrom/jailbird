#!/bin/sh

#run python dose anaysis tool
cellL=A375
timeP=6H

#plate control
#inPut=/xchip/obelix/pod/brew/pc/PRISM001_${cellL}_${timeP}/by_pert_id_pert_dose/PRISM001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n59x978.gctx
inPut=/xchip/obelix/pod/brew/pc/PRISM001_${cellL}_${timeP}/by_pert_id_pert_dose/PRISM001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n60x978.gctx
brewFolder=/xchip/obelix/pod/brew/pc/PRISM001_${cellL}_${timeP}
outDir=/xchip/cogs/hogstrom/analysis/scratch/prism/$cellL/$timeP/ec50_tool

#vehicle control
# inPut=/xchip/obelix/pod/brew/vc/PR_${cellL}_${timeP}/by_pert_id_pert_dose/ASG001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n85x978.gctx
# outDir=/xchip/cogs/projects/ASG_dose_time/cmap_queries/$cellL/$timeP/${metric}_vc

mkdir -p $outDir

python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/ec50_tool.py -r $inPut -o $outDir

#create a broad html page from the results
# ln -s /xchip/cogs/hogstrom/analysis/scratch/ASG_qqPlots/ 

# ln -s /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/jan14/dose_analysis_tool.1358189363930/ PC3_6H_dose