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
	return parser

def main(args,args_from_ui=None):
	'''
	the main entry point into dose_analysis_tool
	'''
	# register the tool with tool ops
	work_dir = tool_ops.register_tool('dose_analysis_tool', args,
	                                   start_path = args.out,
	                                   args_from_ui = args_from_ui)
	

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
	
	
	
	
	###
	## make an SC object from the given gctx file
	#sco = sc.SC()
	#sco.add_sc_from_gctx_meta(args.res, verbose=False)
	#sco.set_thresh_by_specificity(0.8)

	## find all of the unique pert_ids in the data
	#perts = [x.split('::')[0].split(':')[1] for x in sco.pid]
	#unique_perts = set(perts)
	#unique_perts.difference_update(set(['DMSO']))

	## grab the dose information
	#dose = [float(x.split('::')[0].split(':')[2]) for x in sco.pid]

	## grab pert_descs
	#desc = [x.split('::')[1].split(':::')[0] for x in sco.pid]

	## write sc plots to file
	#for unique_pert in unique_perts:
		#sco.plot(include=unique_pert,size=dose,title=unique_pert,
				 #out=os.path.join(work_dir,'_'.join([unique_pert,'SC.png'])))

	## write SC summary table
	#with open(os.path.join(work_dir,'SC_summary.txt'),'w') as f:
		#headers = ['pert_id','pert_desc','base_dose','base_ss',
				   #'base_cc','best_dose','best_ss','best_cc',
				   #'best_ss_lfc','best_cc_lfc','best_sc_lfc_distance']
		#f.write('\t'.join(headers) + '\n')
		#for unique_pert in unique_perts:
			#pert_inds = [i for i,x in enumerate(perts) if unique_pert in x]
			#pert_dose = [dose[x] for x in pert_inds]
			#pert_desc = desc[pert_inds[0]]
			#pert_ss = [sco.s[x] for x in pert_inds]
			#pert_cc = [sco.c[x] for x in pert_inds]
			
			#base_dose = numpy.min(pert_dose)
			#base_ind = pert_dose.index(base_dose)
			#base_ss = pert_ss[base_ind]
			#base_cc = pert_cc[base_ind]
			
			#ss_ratio = numpy.log(numpy.array(pert_ss)/base_ss)
			#cc_ratio = numpy.log(numpy.array(pert_cc)/base_cc)
			#sc_distance = (ss_ratio**2 + cc_ratio**2)**.5
			#sc_distance = sc_distance.tolist()
			
			#best_ind = sc_distance.index(numpy.max(sc_distance))
			#best_dose = pert_dose[best_ind]
			#best_ss = pert_ss[best_ind]
			#best_cc = pert_cc[best_ind]
			#best_ss_ratio = ss_ratio[best_ind]
			#best_cc_ratio = cc_ratio[best_ind]
			#best_sc_distance = sc_distance[best_ind]

			#data = [unique_pert,pert_desc,str(base_dose),str(base_ss),
					#str(base_cc),str(best_dose),str(best_ss),str(best_cc),
					#str(best_ss_ratio),str(best_cc_ratio),str(best_sc_distance)]
			#f.write('\t'.join(data) + '\n')

	## if required, write single probe plots to file
	#if args.probe:
		#gcto = gct.GCT()
		#probe_ind = gcto.get_gctx_rid_inds(args.res,match_list=args.probe,exact=True)
		#gcto.read_gctx_matrix(args.res,row_inds=probe_ind)
		#cids = gcto.get_gctx_cid(args.res)
		#doses = [float(x.split(':')[2]) for x in cids]
		#CM = mu.CMapMongo()
		#with open(os.path.join(work_dir,args.probe + '_summary.txt'),'w') as f:
			#headers = ['pert_id','pert_desc','base_dose','base_z_score',
					   #'best_dose','best_z_score', 'best_z_score_delta']
			#f.write('\t'.join(headers) + '\n')
			#for unique_pert in unique_perts:
				#cid_inds = [i for i,x in enumerate(cids) if unique_pert in x]
				#pert_scores = gcto.matrix[0,cid_inds]
				#pert_doses = [doses[x] for x in cid_inds]
				#tmp_tup = zip(pert_doses,pert_scores)
				#tmp_tup.sort()
				#pert_doses,pert_scores = zip(*tmp_tup)
				#plt.plot(pert_doses,pert_scores)
				#plt.title('::'.join([unique_pert,args.probe]))
				#plt.xlabel('dose')
				#plt.ylabel('z-score')
				#plt.savefig(os.path.join(work_dir,'_'.join([unique_pert,args.probe,'dose_curve.png'])))
				#plt.close()
				
				#pert_desc = CM.find({'pert_id':unique_pert},{'pert_desc':True},limit=1)
				#if not pert_desc:
					#pert_desc = ['-666']
				#pert_desc = pert_desc[0]

				#base_dose = pert_doses[0]
				#base_z_score = pert_scores[0]

				#z_delta = (numpy.array(pert_scores) + 10) - (base_z_score + 10)
				#abs_z_delta = numpy.abs(z_delta)
				#z_delta =  z_delta.tolist()
				#abs_z_delta = abs_z_delta.tolist()
				
				#best_lfc_ind = abs_z_delta.index(numpy.max(abs_z_delta))
				#best_dose = pert_doses[best_lfc_ind]
				#best_z_score = pert_scores[best_lfc_ind]
				#best_z_score_lfc = z_delta[best_lfc_ind]

				#data = [unique_pert,pert_desc,str(base_dose),str(base_z_score),
						#str(best_dose),str(best_z_score),str(best_z_score_lfc)]
				#f.write('\t'.join(data) + '\n')


if __name__ == '__main__':
    # parse the command line arguments
    
    parser = build_parser()
    args = parser.parse_args()
    main(args)