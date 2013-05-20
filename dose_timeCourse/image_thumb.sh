#!/bin/sh

indir=/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/matlab_plots_24H_vc
#indir=/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/ssCCplots_6H
outdir=$indir/thumb600
mkdir $outdir

imgExt=jpg
cd $indir
#loop through the folder and create thumbnales for each image
for i in $(ls *.${imgExt}); do
	convert -size 600x600 $i -resize 600x600 $outdir/$i
done
