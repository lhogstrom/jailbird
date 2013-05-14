// run from the command line as:
// mongo -u cmap_user -p l1000 vitalstatistix/affogato dos_query.js > DOS_sig_ids.txt

// limit query to only DOS compounds, is gold, trt_cp
db.signature.find({sig_id:/DOS/,is_gold:1.0,pert_type:'trt_cp'},{sig_id:1}).forEach(function(x) {print(x.sig_id)});

