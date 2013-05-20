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
	parser.add_argument('-p','--probe', type = str,
	                    default = None, 
	                    help = 'an optional probe to analyze individually')
	parser.add_argument('-o','--out', type = str,
	                    default = os.getcwd(), 
	                    help = 'The output file path')
	parser.add_argument('-q','--query', type = str,
	                    default = 'eslm', 
	                    help = 'metric to be used in cmap query_tool: eslm (enrichment score in landmark space), wteslm (weighted enrichment score), es (enrichment score 22k space)')
	parser.add_argument('-r','--result', action="store_true",
	                    help = 'use this flag to analyze query_tool output after running the query')
	parser.add_argument('-rd','--resultDir', type = str,
	                    help = 'Point dose_tool to a pre-existing directory from query_tool')
	#Rsltgroup = parser.add_mutually_exclusive_group() #two mutually exlusive options: analyze results from 1) working directory or 2) specified directory
	#Rsltgroup.add_argument('-r','--result', action="store_true",
	                    #help = 'use this flag to analyze query_tool output after running the query')
	#Rsltgroup.add_argument('-rd','--resultDir', type = str,
	                    #default = os.getcwd(),
	                    #help = 'Point dose_tool to a pre-existing directory from query_tool')
	return parser

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
	#try:
		#args.resultDir
	#except NameError:
		#outGctx = glob.glob(os.path.join(work_dir, '*COMBINED*.gctx')) #select combined result gctx in working dir created from build_query step
		##rslt.read(outGctx[0])
		#print 'read gct from working dir'
	#else:
		#print args.resultDir
		#print args.result
		##rslt.read(args.resultDir)
		#print 'read gct from explicitly stated result dir'

	if args.result:	
		outGctx = glob.glob(os.path.join(work_dir, '*COMBINED*.gctx')) #select combined result gctx in working dir created from build_query step
		#rslt.read(outGctx[0])
		print 'read gct from working dir'
		print args.result
	else:
		print args.resultDir
		#print args.result
		#rslt.read(args.resultDir)
		print 'read gct from explicitly stated result dir'



	#rsltSigID = rslt.get_rids() #sig IDs from result file

	#qPert = db.get_column_meta('pert_desc')
	#qPertID = db.get_column_meta('pert_id')
	#qDose = db.get_column_meta('pert_dose')
	#ESmat = rslt.matrix
	#iES = ESmat.argsort(axis=0)[::-1] #sort ascending
	#n_inst = len(iES[:,1])

	##loop through each of the perts - graph ranks of query
	#avRnk = []
	#medRnk = []
	#for i, x in enumerate(qPert):
		#iE = iES[:,i] #ES sort index for one column
		#sSigID = []
		#for y in iE:
			#sSigID.append(rsltSigID[y]) #make sorted sig ID list
		#qStr = qPertID[i]
		#cmpd1 = x
		#dose1 = qDose[i]
		#if len(qStr) >= 13:
			#qStr = qStr[0:13] #shorten qPertID
		##i1 = IDsorted.index(qStr) #give first index of match

		##run pymongo query
		#CM = mu.CMapMongo()
		##cmpdSigIds = CM.find({'pert_id':qStr},{'sig_id':True})
		#cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db

		##i1 = all_indices(qStr,sSigID)
		#i1 = [sSigID.index(y) for y in cmpdSigIds] #where instances of the compound of interest sit on the rank list
		#if len(i1) < 1:
			#print cmpd1 + ' has no instances in the cmap database'
			#continue
		#i2 = numpy.array(i1) #convert list to numpy array
		#avr = sum(i2)/len(i2) #what is the average ES rank
		#md = numpy.median(i2) # what is the median ES rank
		#nAv = float(avr)/n_inst #normalize acording to number of instances in db
		#nMd = float(md)/len(iES[:,1]) #normalized median
		#avRnk.append(nAv) #store average ES rank
		#medRnk.append(nMd)
		##plot
		#fname = cmpd1 + '_' + dose1 + '_query_rank.png'
		#outf = os.path.join(work_dir,fname)
		#fig = plt.figure(figsize=(8.0, 2.0))
		#ax = fig.add_subplot(111)
		## the histogram of the data
		#n, bins, patches = ax.hist(i2, 30, facecolor='green', alpha=0.75)
		##ax.set_xlim(0, n_inst)
		#ax.set_xlim(0, int(round(n_inst,-5))) #round instances to nearest 100k
		#ax.set_xlabel('query rank')
		#ax.set_ylabel('freq')
		#ax.set_title('dose = '+ str(dose1) +'um')
		#ax.grid(True)
		#plt.savefig(outf, bbox_inches=0)


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
	#template_heatmap(args,work_dir)
	if args.probe:
		build_probe_curves_and_summary(args,work_dir)
	if args.result or args.resultDir:
		analyze_query(args,work_dir)


def all_indices(value, qlist):
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