#!/bin/sh

OUTPATH=/xchip/cogs/projects/ASG_dose_time/cmap_queries

#list unique compounds on ASG plates
CMPD_LIST=(valproic-acid \
thioridazine \
trifluoperazine \
fulvestrant \
trichostatin-a \
alpha-estradiol \
wortmannin \
tretinoin \
vorinostat \
genistein \
sirolimus \
geldanamycin \
estradiol \
LY-294002 \
alvespimycin \
withaferin-a \
radicicol \
troglitazone \
tanespimycin \
fluphenazine \
mitoxantrone)

#count the number of models
n=${#CMPD_LIST[@]}
# echo $n
#loop through each of the models
#from 0 to n-1 models
for ((s=0; s<n; s++)); do
	echo ${CMPD_LIST[$s]}
	CMPD=${CMPD_LIST[$s]}
	#MCF7
	cp $OUTPATH/MCF7_sigs/matlab_plots_24H/thumb600/${CMPD}_top_genes.jpg $OUTPATH/reports/MCF7_24H_${CMPD}_top_genes.jpg
	cp $OUTPATH/MCF7_sigs/matlab_plots_24H/thumb600/${CMPD}_ss.jpg $OUTPATH/reports/MCF7_24H_${CMPD}_ss.jpg
	cp $OUTPATH/MCF7_sigs/matlab_plots_6H/thumb600/${CMPD}_ss.jpg $OUTPATH/reports/MCF7_6H_${CMPD}_ss.jpg
	cp $OUTPATH/MCF7_sigs/matlab_plots_6H/thumb600/${CMPD}_top_genes.jpg $OUTPATH/reports/MCF7_6H_${CMPD}_top_genes.jpg

	#PC3
	cp $OUTPATH/PC3_sigs/matlab_plots_24H/thumb600/${CMPD}_top_genes.jpg $OUTPATH/reports/PC3_24H_${CMPD}_top_genes.jpg
	cp $OUTPATH/PC3_sigs/matlab_plots_24H/thumb600/${CMPD}_ss.jpg $OUTPATH/reports/PC3_24H_${CMPD}_ss.jpg
	cp $OUTPATH/PC3_sigs/matlab_plots_6H/thumb600/${CMPD}_ss.jpg $OUTPATH/reports/PC3_6H_${CMPD}_ss.jpg
	cp $OUTPATH/PC3_sigs/matlab_plots_6H/thumb600/${CMPD}_top_genes.jpg $OUTPATH/reports/PC3_6H_${CMPD}_top_genes.jpg

	#PC3
	cp $OUTPATH/PC3_sigs/alt_ssCCplots_24H/thumb600/$CMPD.png $OUTPATH/reports/PC3_24H_${CMPD}_ssCCplot.jpg
	cp $OUTPATH/PC3_sigs/alt_ssCCplots_6H/thumb600/$CMPD.png $OUTPATH/reports/PC3_6H_${CMPD}_ssCCplot.jpg

	#MCF7
	cp $OUTPATH/MCF7_sigs/alt_ssCCplots_24H/thumb600/$CMPD.png $OUTPATH/reports/MCF7_24H_${CMPD}_ssCCplot.png
	cp $OUTPATH/MCF7_sigs/alt_ssCCplots_6H/thumb600/$CMPD.png $OUTPATH/reports/MCF7_6H_${CMPD}_ssCCplot.png
done