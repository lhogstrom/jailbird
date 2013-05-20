#! /usr/bin/env python
'''
dose_analysis_tool.py

inputs a gctx file and looks for all uniques perts in the data.  For each pert, dose
information is catalogued and SC plots are made across doses for each pert.  Optionally,
plots for single genes of interest are also made.  Summary data tables are written to disk
as well.

AUTHOR: Corey Flynn, Broad Institute, 2012
'''

import os
import argparse
import numpy
import scipy.stats as stats
import random
import matplotlib.pyplot as plt
import cmap.analytics.sc as sc
import cmap.io.gct as gct
import cmap.util.tool_ops as tool_ops
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import cmap.analytics.signature_strength as ss
import glob


def build_parser():
	'''
	builds the command line parser to use with this tool
	'''
	parser = argparse.ArgumentParser(description = '''inputs a gctx file and looks for all uniques perts in the data.  For each pert, dose
	information is catalogued and SC plots are made across doses for each pert.  Optionally,
	plots for a single probeset of interest are also made.  Summary data tables are written to disk
	as well.''')
	parser.add_argument('res', type = str, help = 'the gctx plate to run analysis on')
	parser.add_argument('-p','--probe', type = str, default = None, help = 'an optional probe to analyze individually')
	parser.add_argument('-o','--out', type = str, default = os.getcwd(),
	                    help = 'The output file path')
	parser.add_argument('-q','--query', type = str, default = 'eslm',
	                    help = 'metric to be used in cmap query_tool: eslm (enrichment score in landmark space), wteslm (weighted enrichment score), es (enrichment score 22k space)')
	Rsltgroup = parser.add_mutually_exclusive_group() #two mutually exlusive options: analyze results from 1) working directory or 2) specified directory
	Rsltgroup.add_argument('-r','--result', action="store_true",
	                    help = 'use this flag to analyze query_tool output after running the query')
	Rsltgroup.add_argument('-rd','--resultDir', type = str,
	                    default = os.getcwd(),
	                    help = 'Point dose_tool to a pre-existing directory from query_tool')
	return parser

def build_SC(args,work_dir):
	'''
	builds SC plots for the dose analysis
	'''
	# instantiate a progress object
	prog = progress.DeterminateProgressBar('Dose Analysis')

	# make an SC object from the given gctx file
	sco = sc.SC()
	sco.add_sc_from_gctx_meta(args.res, verbose=False)
	sco.set_thresh_by_specificity(0.8)

	# find all of the unique pert_ids in the data
	#perts = [':'.join(x.split('::')[0].split(':')[0:2]) for x in sco.pid] $perts is pert_id
	perts = [x.split(':::')[0].split('::')[1] for x in sco.pid] #perts is pert_desc
	pert_ids = [x.split(':')[1] for x in sco.pid]
	unique_perts = set(perts)
	ctl_perts = []
	for i, unique_pert in enumerate(unique_perts):
		#pert_id = unique_pert.split(':')[1]
		#if pert_id == 'DMSO' or pert_id =='CMAP-000':
			#ctl_perts.append(unique_pert)
		if unique_pert == 'DMSO':
			ctl_perts.append(unique_pert)
	unique_perts.difference_update(set(ctl_perts))

	# grab the dose information
	dose = [float(x.split('::')[0].split(':')[2]) for x in sco.pid]

	# grab pert_descs
	desc = [x.split('::')[1].split(':::')[0] for x in sco.pid]

	# write sc plots to file
	num_perts = len(unique_perts)
	for i,unique_pert in enumerate(unique_perts):
		prog.update('making SC plots',i,num_perts)
		sco.plot(include=unique_pert,size=dose,title=unique_pert,pos_con=['None'],out=os.path.join(work_dir,'_'.join([unique_pert.replace(':','_'),'SC.png'])))

	# write SC summary table
	with open(os.path.join(work_dir,'SC_summary.txt'),'w') as f:
		headers = ['pert_id','pert_desc','base_dose','base_ss',
				   'base_cc','best_dose','best_ss','best_cc',
				   'best_ss_lfc','best_cc_lfc','best_sc_lfc_distance']
		f.write('\t'.join(headers) + '\n')
		for i,unique_pert in enumerate(unique_perts):
			prog.update('making SC summary',i,num_perts)
			pert_inds = [i for i,x in enumerate(perts) if unique_pert in x]
			pert_dose = [dose[x] for x in pert_inds]
			pert_desc = desc[pert_inds[0]]
			pert_ss = [sco.s[x] for x in pert_inds]
			pert_cc = [sco.c[x] for x in pert_inds]
			pert_cc = [x if x != -666 else 0 for x in pert_cc]
			
			base_dose = numpy.min(pert_dose)
			base_ind = pert_dose.index(base_dose)
			base_ss = pert_ss[base_ind]
			base_cc = pert_cc[base_ind]
			
			ss_ratio = numpy.log(numpy.array(pert_ss)/base_ss)
			cc_ratio = numpy.log((numpy.array(pert_cc)+1)/(base_cc +1))
			sc_distance = (ss_ratio**2 + cc_ratio**2)**.5
			sc_distance = sc_distance.tolist()
			
			best_ind = sc_distance.index(numpy.max(sc_distance))
			best_dose = pert_dose[best_ind]
			best_ss = pert_ss[best_ind]
			best_cc = pert_cc[best_ind]
			best_ss_ratio = ss_ratio[best_ind]
			best_cc_ratio = cc_ratio[best_ind]
			best_sc_distance = sc_distance[best_ind]

			data = [unique_pert,pert_desc,str(base_dose),str(base_ss),
					str(base_cc),str(best_dose),str(best_ss),str(best_cc),
					str(best_ss_ratio),str(best_cc_ratio),str(best_sc_distance)]
			f.write('\t'.join(data) + '\n')

def build_probe_curves_and_summary(args,work_dir):
	'''
	builds dose response curves for each for the specified probe
	'''
	# instantiate a progress object
	prog = progress.DeterminateProgressBar('Dose Analysis')

	# read the specified probe from the input gctx file
	gcto = gct.GCT()
	probe_ind = gcto.get_gctx_rid_inds(args.res,match_list=args.probe,exact=True)
	gcto.read_gctx_matrix(args.res,row_inds=probe_ind)

	# grab the cids from the file and mine dose information from them.  Find all of 
	# the unique perts
	cids = gcto.get_gctx_cid(args.res)
	doses = [float(x.split(':')[2]) for x in cids]
	perts = [x.split(':')[1] for x in cids]
	unique_perts = list(set(perts))
	
	# for each unique pert_id, find the dose that deviates from the base dose the most.
	# Do template matching to prototype curves. Output a report
	num_perts = len(unique_perts)
	CM = mutil.CMapMongo()
	with open(os.path.join(work_dir,args.probe + '_summary.txt'),'w') as f:
		headers = ['pert_id','pert_desc','base_dose','base_z_score',
				   'best_dose','best_z_score', 'best_z_score_delta',
				   'linear','log','half-log','quarter-log','called shape']
		f.write('\t'.join(headers) + '\n')
		for i,unique_pert in enumerate(unique_perts):
			prog.update('analyzing {0}'.format(args.probe),i,num_perts)
			
			# grab the z-scores and doses for the current pert and sort the pairs
			# by dose
			cid_inds = [i for i,x in enumerate(cids) if unique_pert in x]
			pert_scores = gcto.matrix[0,cid_inds]
			pert_doses = [doses[x] for x in cid_inds]
			tmp_tup = zip(pert_doses,pert_scores)
			tmp_tup.sort()
			pert_doses,pert_scores = zip(*tmp_tup)

			# build the dose response plot for the current pert and save it to disk
			plt.plot(pert_doses,pert_scores)
			plt.title('::'.join([unique_pert,args.probe]))
			plt.xlabel('dose')
			plt.ylabel('z-score')
			plt.savefig(os.path.join(work_dir,'_'.join([unique_pert.replace(':','_'),args.probe,'dose_curve.png'])))
			plt.close()

			# grab the pert_desc from mongo
			pert_desc = CM.find({'pert_id':unique_pert},{'pert_desc':True},limit=1)
			if not pert_desc:
				pert_desc = ['-666']
			pert_desc = pert_desc[0]

			# find the best dose and cast them to lists
			base_dose = pert_doses[0]
			base_z_score = pert_scores[0]

			z_delta = (numpy.array(pert_scores) + 10) - (base_z_score + 10)
			abs_z_delta = numpy.abs(z_delta)
			z_delta =  z_delta.tolist()
			abs_z_delta = abs_z_delta.tolist()
			
			best_ind = z_delta.index(numpy.min(z_delta))
			best_dose = pert_doses[best_ind]
			best_z_score = pert_scores[best_ind]
			best_z_score_delta = z_delta[best_ind]

			if len(pert_doses) > 1:
				# build prototype curves if there is more than one dose
				linear = numpy.linspace(1,10,len(pert_doses))
				log_gen = _log_gen(1)
				log_curve = [log_gen.next() for x in range(len(pert_doses))]
				log_gen = _log_gen(.5)
				half_log_curve = [log_gen.next() for x in range(len(pert_doses))]
				log_gen = _log_gen(.25)
				quarter_log_curve = [log_gen.next() for x in range(len(pert_doses))]

				curves = numpy.array([linear,log_curve,
									  half_log_curve,quarter_log_curve])

				# get the correlation coeficient for each of the curves and the
				# current pert dose curve
				corrs = numpy.corrcoef(pert_scores,curves)
				linear_corr = corrs[0][1]
				log_corr = corrs[0][2]
				half_log_corr = corrs[0][3]
				quarter_log_corr = corrs[0][4]

				#report the best shape by finding the best absolute correlation
				abs_corr = numpy.abs(corrs[0][1:])
				if numpy.where(abs_corr > .8)[0].size > 0:
					abs_corr_max = max(abs_corr)
					abs_corr_max_ind = numpy.where(abs_corr == abs_corr_max)[0][0]
					curve_names = ['linear','log','half-log','quarter-log']
					max_curve_name = curve_names[abs_corr_max_ind]
				else:
					max_curve_name = 'none'

			else:
				# if there is only one dose, set all corrs to 'nan'
				linear_corr = 'nan'
				log_corr = 'nan'
				half_log_corr = 'nan'
				quarter_log_corr = 'nan'
				max_curve_name = 'none'



			# write the dose data to the summary file
			data = [unique_pert,pert_desc,str(base_dose),str(base_z_score),
					str(best_dose),str(best_z_score),str(best_z_score_delta),
					str(linear_corr),str(log_corr),str(half_log_corr),
					str(quarter_log_corr),max_curve_name]
			f.write('\t'.join(data) + '\n')
	prog.clear()

def template_heatmap(args,work_dir):
	'''
	uses template matching to find the most does responsive probesets for each compound in
	the dataset and generates a list of the top 50 and bottom 50 most dose responsive probes.
	heatmaps across all of the doses are made using these probesets
	'''
	# instantiate a progress object
	prog = progress.DeterminateProgressBar('Template Heatmaps')

	# read the data
	gcto = gct.GCT(args.res)
	gcto.read()

	# grab the cids from the file and mine dose information from them.  Find all of 
	# the unique perts
	cids = gcto.get_gctx_cid(args.res)
	pert_descs = gcto.get_column_meta('pert_desc')
	doses = [float(x.split(':')[2]) for x in cids]
	perts = [x.split(':')[1] for x in cids]
	unique_perts = list(set(perts))

	# grab the rid for use below
	rids = gcto.get_gctx_rid(args.res)

	num_perts = len(unique_perts)
	for i,unique_pert in enumerate(unique_perts):
		prog.update('analyzing {0}'.format(unique_pert),i,num_perts)

		# grab the z-scores and doses for the current pert and sort the pairs
		# by dose. put the cid_inds in the same sorted order
		cid_inds = [i for i,x in enumerate(cids) if unique_pert in x]
		pert_desc = pert_descs[cid_inds[0]] #set pert desc to the first dose
		pert_doses = [doses[x] for x in cid_inds]
		tmp_tup = zip(pert_doses,cid_inds)
		tmp_tup.sort()
		pert_doses,cid_inds = zip(*tmp_tup)

		if len(pert_doses) > 1:
			# build prototype curves if there is more than one dose
			linear = numpy.linspace(1,10,len(pert_doses))
			log_gen = _log_gen(1)
			log_curve = [log_gen.next() for x in range(len(pert_doses))]
			log_gen = _log_gen(.5)
			half_log_curve = [log_gen.next() for x in range(len(pert_doses))]
			log_gen = _log_gen(.25)
			quarter_log_curve = [log_gen.next() for x in range(len(pert_doses))]

			curves = numpy.array([linear,log_curve,
								  half_log_curve,quarter_log_curve])

			# correlate all of the probes in the data to the prototype curves
			pert_data = gcto.matrix[:,cid_inds]
			num_probes = pert_data.shape[0]
			cc = numpy.corrcoef(pert_data,curves)

			# grab the correlation values for all the probes against prototype curves
			linear_probe_corrs = cc[0:num_probes,num_probes]
			log_probe_corrs = cc[0:num_probes,num_probes + 1]
			half_log_probe_corrs = cc[0:num_probes,num_probes + 2]
			quarter_log_probe_corrs = cc[0:num_probes,num_probes + 3]

			# compute the random correlation profile for this pert
			num_probes = gcto.matrix.shape[0]
			probe_inds = range(num_probes)
			linear_perm_cc = []
			log_perm_cc = []
			half_log_perm_cc = []
			quarter_log_perm_cc = []
			for i in range(1000):
				perm_curve_inds = [random.sample(probe_inds,1)[0] for x in range(len(pert_doses))]
				perm_curve = [pert_data[perm_curve_inds[x],x] for x in range(len(pert_doses))]
				perm_covar = numpy.corrcoef(perm_curve,curves)
				linear_perm_cc.append(perm_covar[0][1])
				log_perm_cc.append(perm_covar[0][2])
				half_log_perm_cc.append(perm_covar[0][3])
				quarter_log_perm_cc.append(perm_covar[0][4])

			# compute the nominal p values for all correlation values
			linear_probe_corrs_p = numpy.array([stats.percentileofscore(linear_perm_cc,x) 
									for x in linear_probe_corrs])
			log_probe_corrs_p = numpy.array([stats.percentileofscore(log_perm_cc,x) 
									for x in log_probe_corrs])
			half_log_probe_corrs_p = numpy.array([stats.percentileofscore(half_log_perm_cc,x) 
									for x in half_log_probe_corrs])
			quarter_log_probe_corrs_p = numpy.array([stats.percentileofscore(quarter_log_perm_cc,x) 
									for x in quarter_log_probe_corrs])

			# write the p values and correlations out to file
			with open(os.path.join(work_dir,pert_desc + '_template_match_summary.txt'),'w') as f:
				f.write('\t'.join(['probeset','linear corr', 'linear p','log corr', 'log p',
					'half-log corr', 'half-log p','quarter-log corr', 'quarter-log p']) + '\n')
				for j in range(len(linear_probe_corrs)):
					f.write('\t'.join([rids[j],str(linear_probe_corrs[j]), str(linear_probe_corrs_p[j])
						,str(log_probe_corrs[j]), str(log_probe_corrs_p[j])
						,str(half_log_probe_corrs[j]), str(half_log_probe_corrs_p[j])
						,str(quarter_log_probe_corrs[j]), str(quarter_log_probe_corrs_p[j])]) + '\n')


			# build the linear heatmap
			linear_probe_corrs_sort_ind = numpy.argsort(linear_probe_corrs_p)[::-1]
			top = pert_data[linear_probe_corrs_sort_ind[0:50],:]
			bot = pert_data[linear_probe_corrs_sort_ind[-50:],:]
			combined = numpy.vstack([top,bot])
			combined_row_normalized =  combined + numpy.abs(numpy.array([numpy.min(combined,1)]).T)
			row_sums = combined_row_normalized.sum(axis=1)
			combined_row_normalized =  combined_row_normalized / row_sums[:,numpy.newaxis]
			plt.imshow(combined_row_normalized,interpolation='nearest',cmap='RdBu')
			plt.axis('off')
			plt.savefig(os.path.join(work_dir,pert_desc + '_linear_heatmap.png'))

			# build the log heatmap
			log_probe_corrs_sort_ind = numpy.argsort(log_probe_corrs_p)[::-1]
			top = pert_data[log_probe_corrs_sort_ind[0:50],:]
			bot = pert_data[log_probe_corrs_sort_ind[-50:],:]
			combined = numpy.vstack([top,bot])
			combined_row_normalized =  combined + numpy.abs(numpy.array([numpy.min(combined,1)]).T)
			row_sums = combined_row_normalized.sum(axis=1)
			combined_row_normalized =  combined_row_normalized / row_sums[:,numpy.newaxis]
			plt.imshow(combined_row_normalized,interpolation='nearest',cmap='RdBu')
			plt.axis('off')
			plt.savefig(os.path.join(work_dir,pert_desc + '_log_heatmap.png'))

			# build the half log heatmap
			half_log_probe_corrs_sort_ind = numpy.argsort(half_log_probe_corrs_p)[::-1]
			top = pert_data[half_log_probe_corrs_sort_ind[0:50],:]
			bot = pert_data[half_log_probe_corrs_sort_ind[-50:],:]
			combined = numpy.vstack([top,bot])
			combined_row_normalized =  combined + numpy.abs(numpy.array([numpy.min(combined,1)]).T)
			row_sums = combined_row_normalized.sum(axis=1)
			combined_row_normalized =  combined_row_normalized / row_sums[:,numpy.newaxis]
			plt.imshow(combined_row_normalized,interpolation='nearest',cmap='RdBu')
			plt.axis('off')
			plt.savefig(os.path.join(work_dir,pert_desc + '_half_log_heatmap.png'))

			# build the quarter log heatmap
			quarter_log_probe_corrs_sort_ind = numpy.argsort(quarter_log_probe_corrs_p)[::-1]
			top = pert_data[quarter_log_probe_corrs_sort_ind[0:50],:]
			bot = pert_data[quarter_log_probe_corrs_sort_ind[-50:],:]
			combined = numpy.vstack([top,bot])
			combined_row_normalized =  combined + numpy.abs(numpy.array([numpy.min(combined,1)]).T)
			row_sums = combined_row_normalized.sum(axis=1)
			combined_row_normalized =  combined_row_normalized / row_sums[:,numpy.newaxis]
			plt.imshow(combined_row_normalized,interpolation='nearest',cmap='RdBu')
			plt.axis('off')
			plt.savefig(os.path.join(work_dir,pert_desc + '_quarter_log_heatmap.png'))

			# clear that progress bar
			prog.clear()


def build_query(args,work_dir):
	'''
	build query results
	'''
	#make signature for each dose
	fup = os.path.join(work_dir,'up_list.gmt')
	fdn = os.path.join(work_dir,'dn_list.gmt')
	open(fup,'w') #overwrite existing grp file
	open(fdn, 'w') #overwrite existing grp file
	n_edge = 50
	db = gct.GCT()
	#db.read(gctfile)
	db.read(args.res)
	cids = db.get_cids()
	pertIDs = [x.split(':')[1] for x in cids]
	doses = [float(x.split(':')[2]) for x in cids]
	perts = db.get_column_meta('pert_desc')
	probes = db.get_rids()
	cellLs = db.get_column_meta('cell_id')
	timePs = db.get_column_meta('pert_time')	
	mtrx = db.matrix #matrix of data from gct file
	#loop through each column of data
	for i,pertID in enumerate(pertIDs):
		profile = mtrx[:,i]
		n_prof = len(profile)
		iprofile = profile.argsort() #indices that sort array
		iprofile = iprofile[::-1] #switch indicies to decend
		sprofile = profile[iprofile]
		itop = iprofile[0:(n_edge)]
		ibot = iprofile[-n_edge:n_prof]
		col_name = perts[i] + '_' + str(doses[i]) + 'um_' + cellLs[i] + '_' + timePs[i]
		ptop = [] 
		pbot = []
		for j,it in enumerate(itop):
			ptop.append(probes[it]) #make probe id list
		for j,ip in enumerate(ibot):
			pbot.append(probes[ip]) #make probe id list
		#write to gmt list 
		with open(fup,'a') as f:
			f.write(col_name + '\t' + col_name + '\t')
			for pt in ptop:
				f.write(pt + '\t')
			f.write('\n')
		with open(fdn,'a') as f:
			f.write(col_name + '\t' + col_name + '\t')
			for pb in pbot:
				f.write(pb + '\t')
			f.write('\n')
	#python system call
	os.chdir(work_dir)
	#cmd = 'rum -q local query_tool --uptag ' + fup + ' --dntag ' + fdn + ' --metric eslm'
	cmd = 'rum -q local query_tool --uptag ' + fup + ' --dntag ' + fdn + ' --metric wteslm --mkdir false'
	os.system(cmd)

def analyze_query(args,work_dir):
	'''
	Analyze the output from query_tool - find self-connections and create graphs
	'''
	#make a gct object
	db = gct.GCT()
	db.read(args.res)

	##load query result - gctx file
	rslt = gct.GCT()
	#if specific result directory is specified, use that - otherwise get gctx from working dir
	if args.result:
		outGctx = glob.glob(os.path.join(work_dir, '*COMBINED*.gctx')) #select combined result gctx in working dir created from build_query step
		rslt.read(outGctx[0])
	else:
		rslt.read(args.resultDir)

	rsltSigID = rslt.get_rids() #sig IDs from result file

	qPert = db.get_column_meta('pert_desc')
	qPertID = db.get_column_meta('pert_id')
	qDose = db.get_column_meta('pert_dose')
	ESmat = rslt.matrix
	iES = ESmat.argsort(axis=0)[::-1] #sort ascending
	n_inst = len(iES[:,1])

	#loop through each of the perts - graph ranks of query
	prog1 = progress.DeterminateProgressBar('creating self-connection graphs')
	avRnk = []
	medRnk = []
	for i, x in enumerate(qPert):
		prog1.update('graphing {0}',i,len(qPert))
		iE = iES[:,i] #ES sort index for one column
		sSigID = []
		for y in iE:
			sSigID.append(rsltSigID[y]) #make sorted sig ID list
		qStr = qPertID[i]
		cmpd1 = x
		dose1 = qDose[i]
		if len(qStr) >= 13:
			qStr = qStr[0:13] #shorten qPertID
		#i1 = IDsorted.index(qStr) #give first index of match

		#run pymongo query
		CM = mutil.CMapMongo()
		#cmpdSigIds = CM.find({'pert_id':qStr},{'sig_id':True})
		cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db

		#i1 = _all_indices(qStr,sSigID)
		i1 = [sSigID.index(y) for y in cmpdSigIds] #where instances of the compound of interest sit on the rank list
		if len(i1) < 1:
			print cmpd1 + ' has no instances in the cmap database'
			continue
		i2 = numpy.array(i1) #convert list to numpy array
		avr = sum(i2)/len(i2) #what is the average ES rank
		md = numpy.median(i2) # what is the median ES rank
		nAv = float(avr)/n_inst #normalize acording to number of instances in db
		nMd = float(md)/len(iES[:,1]) #normalized median
		avRnk.append(nAv) #store average ES rank
		medRnk.append(nMd)
		#plot
		fname = cmpd1 + '_' + dose1 + '_query_rank.png'
		outf = os.path.join(work_dir,fname)
		fig = plt.figure(figsize=(8.0, 2.0))
		ax = fig.add_subplot(111)
		# the histogram of the data
		n, bins, patches = ax.hist(i2, 30, facecolor='green', alpha=0.75)
		#ax.set_xlim(0, n_inst)
		ax.set_xlim(0, int(round(n_inst,-5))) #round instances to nearest 100k
		ax.set_xlabel('query rank')
		ax.set_ylabel('freq')
		ax.set_title('dose = '+ str(dose1) +'um')
		ax.grid(True)
		plt.savefig(outf, bbox_inches=0)

def external_qq(args,work_dir):
	'''
	make a qq plot of each unique instance - plot the size of each probe acording to how
	often it occurs in the affogato top/bottom 50 list
	'''
	#make a gct object
	db = gct.GCT()
	db.read(args.res)

	qPert = db.get_column_meta('pert_desc')
	qPertID = db.get_column_meta('pert_id')
	qDose = db.get_column_meta('pert_dose')
	probeIDs = db.get_row_meta('id')

	#set null distirbution of z-scores (currently normal)
	ESmat = db.matrix
	#calculate null distribution
	mu, sigma = 0, 1
	s = numpy.random.normal(mu, sigma, len(ESmat[:,1]))
	s.sort()

	pertSet = set(qPert)
	for pert in pertSet:
		iP = _all_indices(pert, qPert) #index of doses on plate
		if len(iP) < 2:
			print pert + ' has only one instance'
			continue
		uDose = [qDose[i] for i in iP]
		fDose = [float(x) for x in uDose] #convert strings to float
		aDose = numpy.asarray(fDose) #convert to numpy array
		iD = aDose.argsort() #local ordering
		sDose = [fDose[j] for j in iD] #sort local doses
		iPo =  [iP[i] for i in iD] #ordered index
		#sMat = ESmat[:,iPo]
		#sMat.sort(axis=0)
		#mongo query for each unique pertID
		qStr = qPertID[iPo[0]] #set pertID
		if len(qStr) >= 13:
			qStr = qStr[0:13] #shorten qPertID
		CM = mutil.CMapMongo()
		#cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True})
		edge50Lst = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True,'up50_lm':True,'dn50_lm':True,'cell_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
		nInstances = len(edge50Lst) #number of instances in db
		#count number of times a probe is in the top/bottom 50 genes of an instance
		upProbeCnts = [0] * len(probeIDs)
		dnProbeCnts = [0] * len(probeIDs)
		for j,inst in enumerate(edge50Lst):
			up50 = edge50Lst[j]['up50_lm']
			dn50 = edge50Lst[j]['dn50_lm']
			#loop through every gene in the top and bottom list - where does it live on the rank list?
			for prb in up50:
				if prb in probeIDs:
					iPrb = probeIDs.index(prb)
					upProbeCnts[iPrb] = upProbeCnts[iPrb] +1
			for prb in dn50:
				if prb in probeIDs:
					iPrb = probeIDs.index(prb)
					dnProbeCnts[iPrb] = dnProbeCnts[iPrb] +1
		#loop through each dose
		for d in iPo:
		#count probe enrichment and plot
				cmpd1 = qPert[d]
				dose1 = qDose[d]
				zLst = db.matrix[:,d]
				iLst = zLst.argsort() #sort z-scores and save index
				sLst = zLst[iLst]
				sUpProbeCnts = [upProbeCnts[l] for l in iLst] #sort probe counts acording to z-score
				sDnProbeCnts = [dnProbeCnts[l] for l in iLst]
				#mkrs = numpy.sqrt(sprobeCnts) # non linear scaling of marker points
				sUpProbeCnts = [float(l) for l in sUpProbeCnts] #convert to float
				sDnProbeCnts = [float(l) for l in sDnProbeCnts] #convert to float
				# upPercMkrs = numpy.divide(sUpProbeCnts,max(sUpProbeCnts)) #divide by max count to make for relative frequency
				# dnPercMkrs = numpy.divide(sDnProbeCnts,max(sDnProbeCnts))
				upPercMkrs = numpy.divide(sUpProbeCnts,nInstances) #divide by total instances to make for relative frequency
				dnPercMkrs = numpy.divide(sDnProbeCnts,nInstances)            
				upMkrs = numpy.multiply(upPercMkrs,100)
				dnMkrs = numpy.multiply(dnPercMkrs,100)
				fig = plt.figure()
				ax = fig.add_subplot(111)
				ax.plot(s,s,'b')
				for j,sl in enumerate(sLst):
					ax.plot(s[j],sl,'r.',markersize=upMkrs[j],alpha=.25)
					ax.plot(s[j],sl,'b.',markersize=dnMkrs[j],alpha=.25)
				ax.set_ylabel('observed z-score')
				ax.set_xlabel('expected z-score')
				# #set legend based on the number of 
				r1 = ax.plot(0,0,'r.',markersize=100,alpha=.25)
				b1 = ax.plot(0,0,'b.',markersize=100,alpha=.25)
				legStrUp = 'probe in 100% of ' + str(nInstances) + ' UP instances'
				legStrDn = 'probe in 100% of ' + str(nInstances) + ' DN instances'
				plt.legend([r1, b1], [legStrUp, legStrDn], numpoints=1, loc=4)
				plt.legdend([b1], ['probe in 100% of ' + str(nInstances) + 'instances' ], numpoints=1)
				ax.set_title(pert + ' dose = ' + dose1)
				fname = pert + '_' + dose1 + 'um_connection_qq.png'
				outf = os.path.join(work_dir,fname)
				plt.savefig(outf, bbox_inches=0)


def main(args,args_from_ui=None):
	'''
	the main entry point into dose_analysis_tool
	'''
	# register the tool with tool ops
	work_dir = tool_ops.register_tool('dose_analysis_tool', args,
	                                   start_path = args.out,
	                                   args_from_ui = args_from_ui)
	# # metics based on gctx input
	# build_SC(args,work_dir)
	# template_heatmap(args,work_dir)
	# if args.probe:
	# 	build_probe_curves_and_summary(args,work_dir)
	# #query based functions
	# if args.query:
	# 	build_query(args,work_dir)
	# if args.result or args.resultDir:
	# 	analyze_query(args,work_dir)
	external_qq(args,work_dir)

def _all_indices(value, qlist):
	'''
	return all indices where the input string matches the item in the list
	'''
	indices = []
	indx = -1
	while True:
		try:
			indx = qlist.index(value, indx+1)
			indices.append(indx)
		except ValueError:
			break
	return indices

def _log_gen(x):
	'''
	generator function to build log series.  x controls the size of the log
	step in the series. x=1 yields full log steps, x=0.5 yields
	half log steps, etc.
	'''
	a = 1
	while True:
		yield a
		a *= 10**x


if __name__ == '__main__':
    # parse the command line arguments
    
    parser = build_parser()
    args = parser.parse_args()
    main(args)