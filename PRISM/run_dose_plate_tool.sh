#!/bin/sh

#run python dose anaysis tool
cellL=MCF7 #PC3, A375, MCF7
timeP=24H

### #plate control
# outDir=/xchip/cogs/projects/PRISM/dose_plate_output-by_pert_id_pert_dose/$cellL/$timeP/$cntrlT
# cntrlT=pc
# inPut=/xchip/cogs/data/brew/a2y13q1/PRISM001_${cellL}_${timeP}/by_pert_id_pert_dose/PRISM001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_*x978.gctx
# brewFolder=/xchip/cogs/data/brew/a2y13q1/PRISM001_${cellL}_${timeP}

#vc control
# cntrlT=vc
# inPut=/xchip/obelix/pod/brew/$cntrlT/PRISM001_${cellL}_${timeP}/by_pert_id_pert_dose/PRISM001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_*x978.gctx
# brewFolder=/xchip/obelix/pod/brew/$cntrlT/PRISM001_${cellL}_${timeP}

#data brewed by_rna_well
cntrlT=pc
inPut=/xchip/cogs/projects/PRISM/data_by_rna_well/pc/PRISM001_${cellL}_${timeP}/by_rna_well/PRISM001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n*x978.gctx
brewFolder=/xchip/cogs/data/brew/a2y13q1/PRISM001_${cellL}_${timeP}
outDir=/xchip/cogs/projects/PRISM/dose_plate_output-by_RNA_well/$cellL/$timeP/$cntrlT



#vehicle control
# inPut=/xchip/obelix/pod/brew/vc/PR_${cellL}_${timeP}/by_pert_id_pert_dose/ASG001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n85x978.gctx
# outDir=/xchip/cogs/projects/ASG_dose_time/cmap_queries/$cellL/$timeP/${metric}_vc

mkdir -p $outDir

python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --brewfolder $brewFolder
# echo /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --brewfolder $brewFolder

#create a broad html page from the results
# ln -s /xchip/cogs/hogstrom/analysis/scratch/ASG_qqPlots/ 

# ln -s /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/jan14/dose_analysis_tool.1358189363930/ /xchip/cogs/web/icmap/hogstrom/PC3_6H_dose