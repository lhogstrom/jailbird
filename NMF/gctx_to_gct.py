#!/bin/sh
use java-1.7

wkdir=/xchip/cogs/projects/NMF/lincs_core_cell_lines


for CELL in A375 A549 HA1E HCC515 HEPG2 HT29 MCF7 PC3 VCAP
do
    file1=$wkdir/$CELL/${CELL}_top_intra_connecting_compound_classes_*.gctx

done