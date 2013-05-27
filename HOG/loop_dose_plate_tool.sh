#!/bin/sh
#run python dose anaysis tool

CELL_LIST=(A549 \
PC3 \
MCF7)
TP_LIST=(6H \
24H)
PLATE_LIST=(HOG001 \
HOG002)

for cellL in $CELL_LIST; do
	for timeP in $TP_LIST; do
		for pName in $PLATE_LIST; do

		echo ${pName}_${cellL}_${timeP}

		done
	done
done

# cellL=A549 #A549, MCF7
# timeP=6H
# 
# #plate control
# cntrlT=pc
# inPut=/xchip/cogs/projects/HOG/data/brew_by_pert_id_pert_dose/pc/HOG001_${cellL}_${timeP}/by_pert_id_pert_dose/HOG001_${cellL}_${timeP}_COMPZ.MODZ_SCORE_LM_n*x978.gctx
# brewFolder=/xchip/cogs/projects/HOG/data/brew_by_pert_id_pert_dose/pc/HOG001_${cellL}_${timeP}
# 
# #vc control
# # cntrlT=vc
# # inPut=
# # brewFolder=/xchip/obelix/pod/brew/$cntrlT/PRISM001_${cellL}_${timeP}
# 
# outDir=/xchip/cogs/projects/HOG/dose_plate_output-by_pert_id_pert_dose/$cellL/$timeP/$cntrlT
# 
# mkdir -p $outDir
# 
# #with brew folder
# # # python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --brewfolder $brewFolder
# # echo /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --brewfolder $brewFolder
# 
# #without brew folder
# python /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --no_summary_self_connect
# # echo /xchip/cogs/hogstrom/scripts/pestle/cmap/tools/dose_plate_tool.py $inPut -o $outDir --no_summary_self_connect
# 
# #create a broad html page from the results
# 
# ln -s $outDir/may24/* /xchip/cogs/web/icmap/hogstrom/HOG_dose_plates/${cellL}_${timeP}