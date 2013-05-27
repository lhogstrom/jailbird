#!/bin/sh
#run python dose anaysis tool
cellL=MCF7 #PC3, A375, MCF7
timeP=24H
cntrlT=pc
platePrefix=HOG003

outDir=/xchip/cogs/projects/HOG/dose_plate_output-by_pert_id_pert_dose/$cellL/$timeP/$cntrlT
rsltDir=/apr25/dose_plate_tool.1366908904964

upFile=${platePrefix}_MCF7_24H_dose_tags_UP_small.gmt
dnFile=${platePrefix}_MCF7_24H_dose_tags_DN_small.gmt

# rum -q local sig_query_tool -metric wtcs -row_space lm -uptag $outDir/$rsltDir/$upFile -dntag $outDir/$rsltDir/$dnFile -out $outDir/$rsltDir
rum -q local sig_query_tool -metric wtcs -row_space lm -uptag $outDir/$rsltDir/$upFile -dntag $outDir/$rsltDir/$dnFile -out $outDir/$rsltDir -column_space gold