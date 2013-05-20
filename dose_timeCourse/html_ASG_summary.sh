#!/bin/sh
###### loop through compound dose responses #######

OUTPATH=/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/dec17/dose_analysis_tool.1355781254466
OUTFILE=cmpd_links.html
#clear contents of outfile
cat /dev/null > $OUTPATH/$OUTFILE

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


#count the number of models
n=${#CMPD_LIST[@]}
# echo $n
#loop through each of the models
#from 0 to n-1 models
for ((s=0; s<n; s++)); do
	echo ${CMPD_LIST[$s]}
	CMPD=${CMPD_LIST[$s]}

	echo "<a href=\"${CMPD}.html\">${CMPD}</a> <BR>"  >> $OUTPATH/$OUTFILE

done