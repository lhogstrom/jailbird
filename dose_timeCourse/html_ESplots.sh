#!/bin/sh
###### loop through compound dose responses #######

OUTPATH=/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_24H_selfConnect_counts/viz
OUTFILE=ASG_ESplots.html
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

	montage ${CMPD}_0.08UM_MCF7_24H_ESscat.jpg  ${CMPD}_0.4UM_MCF7_24H_ESscat.jpg  ${CMPD}_2UM_MCF7_24H_ESscat.jpg  ${CMPD}_10UM_MCF7_24H_ESscat.jpg -tile x4 -geometry 800x200 ${CMPD}_montage_ESscat.jpg

	echo "<h1><CENTER>compound = $CMPD<CENTER></h1> <tr><td><img src=${CMPD}_montage_ESscat.jpg></td></tr>" >> $OUTPATH/$OUTFILE

	montage ${CMPD}_0.08UM_MCF7_24H_EShist.jpg ${CMPD}_0.4UM_MCF7_24H_EShist.jpg  ${CMPD}_2UM_MCF7_24H_EShist.jpg ${CMPD}_10UM_MCF7_24H_EShist.jpg -tile x4 -geometry 800x200 ${CMPD}_montage_EShist.jpg

	echo "<tr><td><img src=${CMPD}_montage_EShist.jpg></td></tr>" >> $OUTPATH/$OUTFILE
done