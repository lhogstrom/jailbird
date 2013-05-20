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
import matplotlib.pyplot as plt
import cmap.analytics.sc as sc
import cmap.io.gct as gct
import cmap.util.tool_ops as tool_ops
import cmap.util.mongo_utils as mu
import cmap.util.progress as progress
import cmap.analytics.signature_strength as ss

def build_parser():
	'''
	builds the command line parser to use with this tool
	'''
	parser = argparse.ArgumentParser(description = '''inputs a gctx file and looks for all uniques perts in the data.  For each pert, dose
	information is catalogued and SC plots are made across doses for each pert.  Optionally,
	plots for a single probeset of interest are also made.  Summary data tables are written to disk
	as well.''')
	parser.add_argument('res', type = str, help = 'the gctx plate to run analysis on')
	parser.add_argument('-p','--probe', type = str,
	                    default = None, 
	                    help = 'an optional probe to analyze individually')
	parser.add_argument('-o','--out', type = str,
	                    default = os.getcwd(), 
	                    help = 'The output file path')
	parser.add_argument('-r','--result', type = str,
	                    default = os.getcwd(), 
	                    help = 'The result directory from query_tool')
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
	unique_perts = set(perts)
	ctl_perts = []
	for unique_pert in unique_perts:
		pert_id = unique_pert.split(':')[1]
		if pert_id == 'DMSO' or pert_id =='CMAP-000':
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

def build_probe_curves(args,work_dir):
	'''
	builds dose response curves for the specified probe
	'''
	gcto = gct.GCT()
	probe_ind = gcto.get_gctx_rid_inds(args.res,match_list=args.probe,exact=True)
	gcto.read_gctx_matrix(args.res,row_inds=probe_ind)
	cids = gcto.get_gctx_cid(args.res)
	doses = [float(x.split(':')[2]) for x in cids]
	CM = mu.CMapMongo()
	with open(os.path.join(work_dir,args.probe + '_summary.txt'),'w') as f:
		headers = ['pert_id','pert_desc','base_dose','base_z_score',
				   'best_dose','best_z_score', 'best_z_score_delta']
		f.write('\t'.join(headers) + '\n')
		for i,unique_pert in enumerate(unique_perts):
			prog.update('analyzing {0}'.format(args.probe),i,num_perts)
			cid_inds = [i for i,x in enumerate(cids) if unique_pert in x]
			pert_scores = gcto.matrix[0,cid_inds]
			pert_doses = [doses[x] for x in cid_inds]
			tmp_tup = zip(pert_doses,pert_scores)
			tmp_tup.sort()
			pert_doses,pert_scores = zip(*tmp_tup)
			plt.plot(pert_doses,pert_scores)
			plt.title('::'.join([unique_pert,args.probe]))
			plt.xlabel('dose')
			plt.ylabel('z-score')
			plt.savefig(os.path.join(work_dir,'_'.join([unique_pert.replace(':','_'),args.probe,'dose_curve.png'])))
			plt.close()
			
			pert_desc = CM.find({'pert_id':unique_pert},{'pert_desc':True},limit=1)
			if not pert_desc:
				pert_desc = ['-666']
			pert_desc = pert_desc[0]

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

			data = [unique_pert,pert_desc,str(base_dose),str(base_z_score),
					str(best_dose),str(best_z_score),str(best_z_score_delta)]
			f.write('\t'.join(data) + '\n')
	prog.clear()

def build_query(args,work_dir):
	'''
	build query results
	'''
	#make signature for each dose
	fup = os.path.join(work_dir,'tmp_up_list.gmt')
	fdn = os.path.join(work_dir,'tmp_dn_list.gmt')
	#fup = '/xchip/cogs/hogstrom/analysis/scratch/tmp_up_list.gmt'
	#fdn = '/xchip/cogs/hogstrom/analysis/scratch/tmp_dn_list.gmt'
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
	cmd = 'rum -q local query_tool --uptag ' + fup + ' --dntag ' + fdn + ' --metric eslm'
	os.system(cmd)

def analyze_query(args,work_dir):
	'''
	Analyze the output from query_tool - make find self-connections and create graphs
	'''
	#make a gct object
	db = gct.GCT()
	db.read(args.res)
	
	#load query result - gctx file
	rslt = gct.GCT()
	rslt.read(args.result)
	
	#import annotations which have been extracted with mongo db
	fname1 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_desc__affogato.txt';
	affPert = open(fname1).read().splitlines()
	fname2 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_id_affogato.txt';
	affPertID = open(fname2).read().splitlines()
	fname3 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/sig_ids_affogato.txt';
	affSigID = open(fname3).read().splitlines()
	rsltSigID = rslt.get_rids()
	
	#sort sigIDs between query result and affogato saved
	iaffSigID = sorted(range(len(affSigID)), key = affSigID.__getitem__) #index of sorting
	irsltSigID = sorted(range(len(rsltSigID)), key = rsltSigID.__getitem__) #index of sorting
	
	#list to be generated according to sigID sort index
	sRid = []
	sAid = []
	sApert = []
	sApertID = []
	for x in irsltSigID:
		sRid.append(rsltSigID[x]) #make sorted sig ID list
	
	for x in iaffSigID:
		sAid.append(affSigID[x]) #make sorted sig ID list
		sApert.append(affPert[x]) #sort perts 
		sApertID.append(affPertID[x]) #sort pertIDs 
	
	if not sAid == sRid:
		print 'result and affogato sig ids are not the same when sorted'
	
	qPert = db.get_column_meta('pert_desc')
	qPertID = db.get_column_meta('pert_id')
	qDose = db.get_column_meta('pert_dose')
	ESmat = rslt.matrix[irsltSigID,:] #order matrix acording to sigID sort
	#iES = ESmat.argsort(axis=0) #sort ascending 
	iES = ESmat.argsort(axis=0)[::-1] #sort ascending 
	n_inst = len(iES[:,1])
	#sortESmat = ESmat[iES]
	
	#loop through each of the perts
	avRnk = []
	medRnk = []
	for i, x in enumerate(qPert):
		iE = iES[:,i] #ES sort index for one column
		IDsorted = []
		for j in iE:
				IDshort = sApertID[j]
				IDshort = IDshort[0:13]
				IDsorted.append(IDshort) 
		qStr = qPertID[i]
		cmpd1 = x
		dose1 = qDose[i]
		if len(qStr) >= 13:
			qStr = qStr[0:13] #shorten qPertID
		#i1 = IDsorted.index(qStr) #give first index of match
		i1 = all_indices(qStr,IDsorted)
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
		ax.set_title('normalized rank score = '+ str(nMd))
		ax.grid(True)
		plt.savefig(outf, bbox_inches=0)


def main(args,args_from_ui=None):
	'''
	the main entry point into dose_analysis_tool
	'''
	# register the tool with tool ops
	work_dir = tool_ops.register_tool('dose_analysis_tool', args,
	                                   start_path = args.out,
	                                   args_from_ui = args_from_ui)
	#build_SC(args,work_dir)
	#build_query(args,work_dir)
	if args.probe:
		build_probe_curves(args,work_dir)	
	if args.result:
		analyze_query(args,work_dir)

#return all indices where the input string matches the item in the list
def all_indices(value, qlist):
	indices = []
	indx = -1
	while True:
		try:
			indx = qlist.index(value, indx+1)
			indices.append(indx)
		except ValueError:
			break
	return indices


if __name__ == '__main__':
    # parse the command line arguments
    
    parser = build_parser()
    args = parser.parse_args()
    main(args)