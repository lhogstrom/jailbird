#!/bin/sh
###### loop through compound dose responses - make unique html page for each one #######

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

	OUTFILE=reports/${CMPD}.html
	#clear contents of outfile
	cat /dev/null > $OUTPATH/$OUTFILE


# 	dose1=awk 'NR == 1 {print; exit}' $OUTPATH/MCF7_6H_selfConnect_counts/thioridazine_query_connections.txt
# 	echo $dose1

	#print python figs	
# 	echo "<h1><CENTER>compound = $CMPD<CENTER></h1> <tr><td><img src=$OUTPATH/PC3_sigs/ssCCplots_6H/thumb600/$CMPD.png></td><td><img src=$OUTPATH/PC3_sigs/ssCCplots_24H/thumb600/$CMPD.png></td></tr><tr><td><img src=$OUTPATH/MCF7_sigs/ssCCplots_6H/thumb600/$CMPD.png></td><td><img src=$OUTPATH/MCF7_sigs/ssCCplots_24H/thumb600/$CMPD.png></td></tr>" >> $OUTPATH/$OUTFILE
	#alt ss calculation	
	echo "<h1><CENTER>compound = $CMPD<CENTER></h1> <tr><td><img src=$OUTPATH/PC3_sigs/alt_ssCCplots_6H/thumb600/$CMPD.png></td><td><img src=$OUTPATH/PC3_sigs/alt_ssCCplots_24H/thumb600/$CMPD.png></td></tr><tr><td><img src=$OUTPATH/MCF7_sigs/alt_ssCCplots_6H/thumb600/$CMPD.png></td><td><img src=$OUTPATH/MCF7_sigs/alt_ssCCplots_24H/thumb600/$CMPD.png></td></tr>" >> $OUTPATH/$OUTFILE
	
	echo "<BR>" >> $OUTPATH/$OUTFILE

	#print matlab figs
	echo "<tr><td><img src=$OUTPATH/PC3_sigs/matlab_plots_6H/thumb600/${CMPD}_top_genes.jpg></td><td><img src=$OUTPATH/PC3_sigs/matlab_plots_24H/thumb600/${CMPD}_top_genes.jpg></td></tr><tr><td><img src=$OUTPATH/MCF7_sigs/matlab_plots_6H/thumb600/${CMPD}_top_genes.jpg></td><td><img src=$OUTPATH/MCF7_sigs/matlab_plots_24H/thumb600/${CMPD}_top_genes.jpg></td></tr>" >> $OUTPATH/$OUTFILE	

	echo "<h2>affogato self-connection events - from top 1000 ES rank</h2>" >> $OUTPATH/$OUTFILE
	#PC3 6H self-connections
 	awk 'BEGIN {print "<table border="1"><tr><th>PC3 6H</th>"} END {print "</table>"}{print "<tr>\n<td>dose "NR" - "$0"</td>\n</tr>"}' /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_6H_selfConnect_counts/${CMPD}_query_connections.txt >> $OUTPATH/$OUTFILE

	#PC3 24H self-connections
 	awk 'BEGIN {print "<table border="1"><tr><th>PC3 24H</th>"} END {print "</table>"}{print "<tr>\n<td>dose "NR" - "$0"</td>\n</tr>"}' /xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_24H_selfConnect_counts/${CMPD}_query_connections.txt >> $OUTPATH/$OUTFILE

	#MCF7 6H self-connections
 	awk 'BEGIN {print "<table border="1"><tr><th>MCF7 6H</th>"} END {print "</table>"}{print "<tr>\n<td>dose "NR" - "$0"</td>\n</tr>"}' /xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_6H_selfConnect_counts/${CMPD}_query_connections.txt >> $OUTPATH/$OUTFILE

	#MCF7 24H self-connections
 	awk 'BEGIN {print "<table border="1"><tr><th>MCF7 24H</th>"} END {print "</table>"}{print "<tr>\n<td>dose "NR" - "$0"</td>\n</tr>"}' /xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_24H_selfConnect_counts/${CMPD}_query_connections.txt >> $OUTPATH/$OUTFILE

done