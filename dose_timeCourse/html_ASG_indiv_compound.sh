#!/bin/sh
###### loop through compound dose responses - make unique html page for each one #######

OUTPATH=/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/dec17/dose_analysis_tool.1355781254466

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

#cherry picked list
# CMPD_LIST=(trichostatin-a \
# wortmannin \
# vorinostat \
# sirolimus \
# geldanamycin \
# LY-294002 \
# alvespimycin \
# withaferin-a \
# radicicol \
# tanespimycin \
# mitoxantrone)

cd $OUTPATH
#count the number of models
n=${#CMPD_LIST[@]}
# echo $n
#loop through each of the models
#from 0 to n-1 models
for ((s=0; s<n; s++)); do
	echo ${CMPD_LIST[$s]}
	CMPD=${CMPD_LIST[$s]}

	OUTFILE=${CMPD}.html
	#clear contents of outfile
	cat /dev/null > $OUTPATH/$OUTFILE

	montage ${CMPD}_0.08_query_rank.png  ${CMPD}_0.40_query_rank.png ${CMPD}_2.00_query_rank.png ${CMPD}_10.00_query_rank.png -tile x4 -geometry 800x200 ${CMPD}_montage_rank.png

	echo "<h1><CENTER>compound = $CMPD<CENTER></h1><h2>affogato query metrics</h2> <tr><td><img src=${CMPD}_montage_rank.png></td></tr>" >> $OUTPATH/$OUTFILE

	echo "<BR>" >> $OUTPATH/$OUTFILE	
	echo "<tr><td><img src=${CMPD}_SC.png></td></tr>" >> $OUTPATH/$OUTFILE
# 	echo "<tr><td><img src=${CMPD}_normalized_rank.png></td></tr>" >> $OUTPATH/$OUTFILE

	echo "<BR>" >> $OUTPATH/$OUTFILE
	echo "<tr><td><img src=${CMPD}_linear_heatmap.png></td></tr>" >> $OUTPATH/$OUTFILE
	echo "<tr><td><img src=${CMPD}_log_heatmap.png></td></tr>" >> $OUTPATH/$OUTFILE
	echo "<tr><td><img src=${CMPD}_half_log_heatmap.png></td></tr>" >> $OUTPATH/$OUTFILE
	echo "<tr><td><img src=${CMPD}_quarter_log_heatmap.png></td></tr>" >> $OUTPATH/$OUTFILE

	echo "<BR>" >> $OUTPATH/$OUTFILE
	echo "<tr><td><img src=${CMPD}_0.08um_internal-external_qq.png></td></tr>" >> $OUTPATH/$OUTFILE
	echo "<tr><td><img src=${CMPD}_0.40um_internal-external_qq.png></td></tr>" >> $OUTPATH/$OUTFILE
	echo "<tr><td><img src=${CMPD}_2.00um_internal-external_qq.png></td></tr>" >> $OUTPATH/$OUTFILE
	echo "<tr><td><img src=${CMPD}_10.00um_internal-external_qq.png></td></tr>" >> $OUTPATH/$OUTFILE

done