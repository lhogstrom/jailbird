#!/bin/sh

OUTPATH=/xchip/cogs/projects/ASG_dose_time/cmap_queries/cmpd_sigIDs

CMPD_LIST=(trifluoperazine)
ID_LIST=(BRD-K89732114)

# CMPD_LIST=(alvespimycin \
# troglitazone \
# trichostatin-a \
# geldanamycin \
# radicicol \
# alpha-estradiol \
# tanespimycin \
# wortmannin \
# sirolimus \
# thioridazine \
# fulvestrant \
# estradiol \
# mitoxantrone \
# LY-294002 \
# valproic-acid \
# genistein \
# fluphenazine \
# tretinoin \
# vorinostat \
# withaferin-a \
# trifluoperazine)
# 
# ID_LIST=(BRD-A06304526 \
# BRD-A13084692 \
# BRD-A19037878 \
# BRD-A19500257 \
# BRD-A39996500 \
# BRD-A60070924 \
# BRD-A61304759 \
# BRD-A75409952 \
# BRD-A79768653 \
# BRD-A84481105 \
# BRD-A90490067 \
# BRD-K18910433 \
# BRD-K21680192 \
# BRD-K27305650 \
# BRD-K41260949 \
# BRD-K43797669 \
# BRD-K55127134 \
# BRD-K71879491 \
# BRD-K81418486 \
# BRD-K88378636 \
# BRD-K89732114)

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
	PERTID=${ID_LIST[$s]}
	
	## search by pert_desc ##
	#echo "db.signature.find({pert_desc:/${CMPD}/,is_gold:1.0},{sig_id:1}).forEach(function(x) {print(x.sig_id)});" >> $OUTPATH/tmp_quer.js

	## search by pert_id - save sig ID ##
	echo "db.signature.find({pert_id:/${PERTID}/},{sig_id:1}).forEach(function(x) {print(x.sig_id)});" > $OUTPATH/tmp_quer.js
	## search by pert_id - save pert ID ##
	#echo "db.signature.find({pert_id:/${PERTID}/},{pert_id:1}).forEach(function(x) {print(x.sig_id)});" > $OUTPATH/tmp_quer.js

	mongo --quiet -u cmap_user -p l1000 vitalstatistix/affogato $OUTPATH/tmp_quer.js > $OUTPATH/${CMPD}_bypertID.grp

done