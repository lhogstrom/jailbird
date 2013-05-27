#!/bin/sh

#run python dose anaysis tool
cellL=MCF7 #PC3, A375, MCF7
timeP=6H

#data brewed by_pert_id_pert_dose
cntrlT=pc
inPut=/xchip/cogs/data/brew/a2y13q1/HDAC002_${cellL}_${timeP}/by_pert_id_pert_dose/HDAC002_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n*x978.gctx
brewFolder=/xchip/cogs/data/brew/a2y13q1/HDAC002_${cellL}_${timeP}
outDir=/xchip/cogs/hogstrom/analysis/HDAC/dose_analysis/${cellL}_${timeP}/$cntrlT

mkdir -p $outDir

# python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --brewfolder $brewFolder
python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --brewfolder $brewFolder --no_SC
# echo /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --brewfolder $brewFolder --no_SC

#create a broad html page from the results
# ln -s /xchip/cogs/hogstrom/analysis/HDAC/dose_analysis/PC3_24H/pc/may24/dose_plate_tool.1369409043173/ /xchip/cogs/web/icmap/hogstrom/HDAC002_dose_analysis/PC3_24H

# ln -s /xchip/cogs/hogstrom/analysis/HDAC/dose_analysis/MCF7_6H/pc/may24/dose_plate_tool.1369409129349/ /xchip/cogs/web/icmap/hogstrom/HDAC002_dose_analysis/MCF7_6H

# ls $outDir/may24/*