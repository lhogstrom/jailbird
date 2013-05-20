// run from the command line as:
// mongo --quiet -u cmap_user -p l1000 vitalstatistix/affogato /xchip/cogs/hogstrom/scripts/dose_timeCourse/compound_query.js > /xchip/cogs/projects/ASG_dose_time/cmap_queries/cmpd_sigIDs/wortmannin.grp

// limit query to only DOS compounds, is gold, trt_cp
// db.signature.find({is_gold:1.0},{sig_id:1}).forEach(function(x) {print(x.sig_id)});

// limit query to only DOS compounds, is gold, trt_cp
// db.signature.find({sig_id:/DOS/,is_gold:1.0,pert_type:'trt_cp'},{sig_id:1}).forEach(function(x) {print(x.sig_id)});

	// search all for all instances of a sepcific compound - by pert_desc
	// db.signature.find({pert_desc:/wortmannin/,is_gold:1.0,pert_type:'trt_cp'},{sig_id:1}).forEach(function(x) {print(x.sig_id)});

	// search all for all instances of a sepcific compound - by pert_id
	db.signature.find({pert_id:/BRD-K87343924-001-02-4/,is_gold:1.0},{sig_id:1}).forEach(function(x) {print(x.sig_id)});