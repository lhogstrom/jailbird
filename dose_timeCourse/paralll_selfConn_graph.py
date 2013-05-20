import os
import jinja2
import argparse
import numpy
import scipy.stats as stats
import random
import matplotlib.pyplot as plt
import cmap
import cmap.analytics.sc as sc
import cmap.io.gct as gct
import cmap.util.tool_ops as tool_ops
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import cmap.analytics.signature_strength as ss
import glob
import subprocess
import time
import re
import operator
import multiprocessing

def _self_connection_graph(params):
	'''
	input: the index of a signature from the input gctx, output: saves a self connection graph
	'''
	cmpd1 = params[0]
	dose1 = params[1]
	i2 = params[2] 
	work_dir = params[3]
	n_inst = params[4]

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
	# indices = [i for i, x in enumerate(qlist) if x == value]
	return indices

# # def analyze_query(args,work_dir):
# 	'''
# 	Analyze the output from query_tool - find self-connections and create graphs
# 	'''
# 	#make a gct object
# 	# db = gct.GCT()
# 	# db.read(args.res)

# work_dir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/scratch'
# fname = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
# db = gct.GCT()
# db.read(fname)

# ##load query result - gctx file
# rslt = gct.GCT()
# #if specific result directory is specified, use that - otherwise get gctx from working dir
# rslt.read('/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/jan04/dose_analysis_tool.1357311906761/result_WTESLM.COMBINED_n85x398050.gctx')

qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qDose = db.get_column_meta('pert_dose')
ESmat = rslt.matrix
iES = ESmat.argsort(axis=0)[::-1] #sort ascending
n_inst = len(iES[:,1])

#loop through each of the perts - graph ranks of query
prog1 = progress.DeterminateProgressBar('creating self-connection graphs')
qSumF = os.path.join(work_dir, 'percent_top_query_self_connections.txt')
avRnk = []
medRnk = []
prRnk = []
tupList = []
#loop through each of the UNIQUE perts - graph ranks of query
pertSet = set(qPert)
for pert in pertSet:
	cmpd1 = pert
	iP = _all_indices(pert, qPert) #index of doses on plate
	if len(iP) < 2:
		print pert + ' has only one instance on the plate'
		continue
	uDose = [qDose[i] for i in iP]
	fDose = [float(x) for x in uDose] #convert strings to float
	aDose = numpy.asarray(fDose) #convert to numpy array
	iD = aDose.argsort() #local ordering
	sDose = [fDose[j] for j in iD] #sort local doses
	iPo =  [iP[i] for i in iD] #ordered index
	qStr = qPertID[iPo[0]] #set pertID
	if len(qStr) >= 13:
		qStr = qStr[0:13] #shorten qPertID
	#run pymongo query
	CM = mutil.CMapMongo()
	#cmpdSigIds = CM.find({'pert_id':qStr},{'sig_id':True})
	cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
	n_cmpdInst = len(cmpdSigIds)
	if len(cmpdSigIds) < 1:
		print cmpd1 + ' has no instances in the cmap database'
		continue
	if len(cmpdSigIds) == 1:
		print cmpd1 + ' has one instance in the cmap database'
	#loop through each dose ()
	topTable = []
	for d in iPo:
	#count probe enrichment and plot
			cmpd1 = qPert[d]
			dose1 = qDose[d]
			iE = iES[:,d] #ES sort index for one column
			sSigID = []
			for y in iE:
				sSigID.append(rsltSigID[y]) #make sorted sig ID list
			i1 = [sSigID.index(y) for y in cmpdSigIds] #where instances of the compound of interest sit on the rank list
			i2 = numpy.array(i1) #convert list to numpy array
			avr = sum(i2)/len(i2) #what is the average ES rank
			md = numpy.median(i2) # what is the median ES rank
			nAv = float(avr)/n_inst #normalize acording to number of instances in db
			nMd = float(md)/len(iES[:,1]) #normalized median
			i1.sort()
			# np = 1000 #define the top of the rank
			np = int(n_inst*.01) # 1% of all instances
			ntop = [x for x in i1 if x <= np]
			# PercTop = ntop/n_cmpdInst #percent of instances in the top of the rank
			nPr = float(len(ntop))/(len(i1)) #percent of instances at the top of the list
			topTable.append(nPr)
			prRnk.append(nPr)
			avRnk.append(nAv) #store average ES rank
			medRnk.append(nMd)
			# make command tuple
			tup1 = (cmpd1,dose1,i2,work_dir,n_inst) #list of tuples
			tupList.append(tup1) #add to list of tuples 
	#top list for each compound 
	with open(qSumF,'a') as f:
		f.write(cmpd1 + '\t')
		for pr in topTable:
			f.write(str(pr) + '\t')
		f.write('\n')

# instantiate a progress object
prog = progress.DeterminateProgressBar('self-connection graph builder')
#build graphs in parallel
pool = multiprocessing.Pool()
rs = pool.map_async(_self_connection_graph,tupList)
pool.close() # No more work
while (True):
    if (rs.ready()): break
    remaining = rs._number_left
    prog.show_message('Waiting for {0} tasks to complete...'.format(remaining))
    time.sleep(0.1)


