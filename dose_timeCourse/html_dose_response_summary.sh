#!/bin/sh
###### loop through compound dose responses #######

# OUTPATH=/xchip/cogs/hogstrom/analysis/scratch/Nov27
OUTPATH=/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7/24H/wteslm/dec10/dose_analysis_tool.1355156388923
# cellLine=MCF7
# timeP=24H
# FPATH2=/xchip/cogs/projects/ASG_dose_time/cmap_queries/reports
OUTFILE=ASG_rank_sum.html
#clear contents of outfile
cat /dev/null > $OUTPATH/$OUTFILE

cd $OUTPATH

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

	montage ${CMPD}_0.08_query_rank.png  ${CMPD}_0.40_query_rank.png ${CMPD}_2.00_query_rank.png ${CMPD}_10.00_query_rank.png -tile x4 -geometry 800x200 ${CMPD}_montage_rank.png

	echo "<h1><CENTER>compound = $CMPD<CENTER></h1><h2>affogato query metrics</h2> <tr><td><img src=${CMPD}_montage_rank.png></td></tr>" >> $OUTPATH/$OUTFILE

	echo "<BR>" >> $OUTPATH/$OUTFILE	
	echo "<tr><td><img src=${CMPD}_SC.png></td></tr>" >> $OUTPATH/$OUTFILE
# 	echo "<tr><td><img src=${CMPD}_normalized_rank.png></td></tr>" >> $OUTPATH/$OUTFILE

	echo "<tr><td><img src=${CMPD}-heatmap.png></td></tr>" >> $OUTPATH/$OUTFILE

# 	echo "<h2>profile signal metrics</h2><tr><td><img src=${CMPD}_ss.png></td></tr>" >> $OUTPATH/$OUTFILE #regular ss
# 
# 	echo "<tr><td><img src=${CMPD}_ss_thresh.png></td></tr>" >> $OUTPATH/$OUTFILE #threshold z > 2 ss

# #output from ss calculations
# 	echo "<BR>" >> $OUTPATH/$OUTFILE
# 	#ss plot
# # 	echo "<tr><td><img src=$FPATH2/MCF7_24H_${CMPD}_ss.jpg></td></tr>" >> $OUTPATH/$OUTFILE
# 	#sc plot
# 	echo "<tr><td><img src=$FPATH2/${cellLine}_${timeP}_${CMPD}_ssCCplot.jpg></td></tr>" >> $OUTPATH/$OUTFILE
# 	#top genes plot
# 	echo "<tr><td><img src=$FPATH2/${cellLine}_${timeP}_${CMPD}_top_genes.jpg></td></tr>" >> $OUTPATH/$OUTFILE
# 	echo "<tr><td><img src=${CMPD}_qq.png></td></tr>" >> $OUTPATH/$OUTFILE #qq plot
done