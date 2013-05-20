#!/bin/sh
###### loop through multiple glm models and create overlay images #######

CTRLDIR=/usit/tsd/nevro-work/projects/DIV/larson/control
DATADIR=/usit/tsd/nevro-work/projects/multisample/5.1/memory
SCRIPTSDIR=/usit/tsd/nevro-work/projects/DIV/larson/scripts

#export freesurfer directory as global variable
export SUBJECTS_DIR=/usit/tsd/nevro-work/projects/multisample/5.1/memory
cd $DATADIR


#define outpath
FWHM=30
THRESH=2
OUTPATH=$DATADIR/model_overaly_images/fwhm$FWHM/age_bin_thresh2_lgi
OUTFILE=20yr_age_bin_model_image_report.html

#### create html header ####
echo "<h2>glm model overlays - fwhm = $FWHM, THRESH = $THRESH</h2>" > $OUTPATH/$OUTFILE



#list second model directory - must match the model in MODEL_LIST1
#MODEL_LIST2=(lh_thickness_pvr_area_regress_out_age_gender)
MODEL_LIST2=(lh_avg_area_cor \
surface_measure_age_regress_totalWMvol_gender \
surface_measure_age_regress_totalWMvol_gender \
surface_measure_w_totalWMvol_regress_age_gender \
surface_measure_w_totalWMvol_regress_age_gender \
lh_thickness_pvr_area_regress_out_age_gender \
lh_thickness_pvr_area_regress_out_age_gender)

	#count the number of models
	n=${#MODEL_LIST2[@]}
	#echo $n
	#loop through each of the models
	#from 0 to n-1 models
	for ((s=0; s<n; s++)); do
	echo ${MODEL_LIST2[$s]}
		
		for bin in young middle old; do 
		AGEBIN=$bin
		
		#list the first model directory (withouth hemisphere)
		#MODEL_LIST1=(age_bin_$AGEBIN.thickness_pvr_white_regress_age_gender.glmdir)
		MODEL_LIST1=(age_bin_$AGEBIN.pial_lgi_w_age_regress_gender.glmdir \
		age_bin_$AGEBIN.pial_lgi_w_age.regress_totSA_gender.glmdir \
		age_bin_$AGEBIN.pial_lgi_w_age.regress_totWMvol_gender.glmdir \
		age_bin_$AGEBIN.pial_lgi_w_totSA.regress_age_gender.glmdir \
		age_bin_$AGEBIN.pial_lgi_w_totWMvol.regress_age_gender.glmdir \
		age_bin_$AGEBIN.thickness_pvr_pial_lgi_regress_age_gender.glmdir \
		age_bin_$AGEBIN.white_pvr_pial_lgi_regress_age_gender.glmdir)

			#loop through both hemispheres
		# 	for hem1 in lh rh; do
			
			#MODEL=$hem1.322.area_white_age.regress_thickness_pvr_gender.glmdir
			MODEL=${MODEL_LIST1[$s]}
		
			#echo "$OUTPATH/rh.$MODEL/fsaverage_wh_rh_lat.png"
		
			echo "<h4>age bin = $AGEBIN - MODEL = $MODEL</h4> <tr><td><img src=$OUTPATH/$AGEBIN/rh.$MODEL/fsaverage_wh_rh_lat_250.png></td><td><img src=$OUTPATH/$AGEBIN/rh.$MODEL/fsaverage_wh_rh_med_250.png></td></tr><tr><td><img src=$OUTPATH/$AGEBIN/lh.$MODEL/fsaverage_wh_lh_med_250.png></td><td><img src=$OUTPATH/$AGEBIN/lh.$MODEL/fsaverage_wh_lh_lat_250.png></td></tr>" >> $OUTPATH/$OUTFILE
		
		# 	done
	
		done

	done