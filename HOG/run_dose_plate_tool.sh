#!/bin/sh
#run python dose anaysis tool
cellL=A549 #A549, MCF7
timeP=6H

### run individual roast plates
# cntrlT=pc
# # inPut=/xchip/obelix/pod/roast/HOG001_${cellL}_${timeP}_X1_B10_DUO52HI53LO/zs/HOG001_${cellL}_${timeP}_X1_B10_DUO52HI53LO_ZSPCQNORM_n374x978.gct
# inPut=/xchip/obelix/pod/roast/HOG002_MCF7_24H_X1_B10_DUO52HI53LO/zs/HOG002_MCF7_24H_X1_B10_DUO52HI53LO_ZSPCQNORM_n369x978.gct

#plate control
cntrlT=pc
inPut=/xchip/cogs/projects/HOG/data/brew_by_pert_id_pert_dose/pc/HOG001_${cellL}_${timeP}/by_pert_id_pert_dose/HOG001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n*x978.gctx
brewFolder=/xchip/cogs/projects/HOG/data/brew_by_pert_id_pert_dose/pc/HOG001_${cellL}_${timeP}

#vc control
# cntrlT=vc
# inPut=
# brewFolder=/xchip/obelix/pod/brew/$cntrlT/PRISM001_${cellL}_${timeP}

#data brewed by_rna_well
# inPut=/xchip/cogs/projects/PRISM/data_by_rna_well/pc/PRISM001_${cellL}_${timeP}/by_rna_well/PRISM001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n260x978.gctx

outDir=/xchip/cogs/projects/HOG/dose_plate_output-by_pert_id_pert_dose/$cellL/$timeP/$cntrlT

mkdir -p $outDir

#with brew folder
# # python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --brewfolder $brewFolder
# echo /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --brewfolder $brewFolder

#without brew folder
python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --no_summary_self_connect
# echo /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --no_summary_self_connect

#create a broad html page from the results

ln -s $outDir/may24/* /xchip/cogs/web/icmap/hogstrom/HOG_dose_plates/${cellL}_${timeP}