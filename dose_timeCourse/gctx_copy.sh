#!/bin/bash

INDIR=/xchip/obelix/pod/roast
OUTDIR=/xchip/cogs/projects/ASG_dose_time/data

cd $INDIR
#ls ASG*/ -d1 >> $OUTDIR/plate_list.txt

#grab QNORM data
# 	cat $OUTDIR/plate_list.txt | while read line; do #loop through each line of the text file
# 		cp ${line}/${line}_INF* $OUTDIR
# 		cp ${line}/${line}_QNORM* $OUTDIR
# 	done

#grab zscore data
	cat $OUTDIR/plate_list.txt | while read line; do #loop through each line of the text file
		cp ${line}/zs/${line}_ZSPCQNORM* $OUTDIR
		cp ${line}/zs/${line}_ZSVCQNORM* $OUTDIR
	done