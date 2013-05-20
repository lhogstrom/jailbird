#/!/bin/py
#run file like this
#execfile('/xchip/cogs/hogstrom/scripts/dose_timeCourse/py_profiler_dose.py')
#run from command line like this: python -m cProfile -o py_profiler.profile py_profiler_dose.py


import cProfile
import os
import jinja2
import argparse
import numpy
import scipy.stats as stats
import random
import matplotlib.pyplot as plt
import cmap
#import cmap.analytics.sc as sc
import cmap.io.gct as gct
import cmap.util.tool_ops as tool_ops
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import cmap.analytics.signature_strength as ss
import glob
import subprocess


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

work_dir = '/xchip/cogs/hogstrom/analysis/scratch/ASG_qqPlots/test1'
cellLine = 'MCF7'
timeP = '24H'
gctfile = '/xchip/obelix/pod/brew/pc/ASG001_%s_%s/by_pert_id_pert_dose/ASG001_%s_%s_COMPZ.MODZ_SCORE_LM_n85x978.gctx' % (cellLine,timeP,cellLine,timeP) 
#make a gct object
db = gct.GCT()
db.read(gctfile)

### Build profiler ####

# #make signature for each dose
# fup = os.path.join(work_dir,'up_list.gmt')
# fdn = os.path.join(work_dir,'dn_list.gmt')
# open(fup,'w') #overwrite existing grp file
# open(fdn, 'w') #overwrite existing grp file
# n_edge = 50
# cids = db.get_cids()
# pertIDs = [x.split(':')[1] for x in cids]
# doses = [float(x.split(':')[2]) for x in cids]
# perts = db.get_column_meta('pert_desc')
# probes = db.get_rids()
# cellLs = db.get_column_meta('cell_id')
# timePs = db.get_column_meta('pert_time')	
# mtrx = db.matrix #matrix of data from gct file
# #loop through each column of data
# for i,pertID in enumerate(pertIDs):
# 	profile = mtrx[:,i]
# 	n_prof = len(profile)
# 	iprofile = profile.argsort() #indices that sort array
# 	iprofile = iprofile[::-1] #switch indicies to decend
# 	sprofile = profile[iprofile]
# 	itop = iprofile[0:(n_edge)]
# 	ibot = iprofile[-n_edge:n_prof]
# 	col_name = perts[i] + '_' + str(doses[i]) + 'um_' + cellLs[i] + '_' + timePs[i]
# 	ptop = [] 
# 	pbot = []
# 	for j,it in enumerate(itop):
# 		ptop.append(probes[it]) #make probe id list
# 	for j,ip in enumerate(ibot):
# 		pbot.append(probes[ip]) #make probe id list
# 	#write to gmt list 
# 	with open(fup,'a') as f:
# 		f.write(col_name + '\t' + col_name + '\t')
# 		for pt in ptop:
# 			f.write(pt + '\t')
# 		f.write('\n')
# 	with open(fdn,'a') as f:
# 		f.write(col_name + '\t' + col_name + '\t')
# 		for pb in pbot:
# 			f.write(pb + '\t')
# 		f.write('\n')
# ## systems call to run query tool ## 
# os.chdir(work_dir)
# # #cmd = 'rum -q local query_tool --uptag ' + fup + ' --dntag ' + fdn + ' --metric eslm'
# # cmd = 'rum -q local query_tool --uptag ' + fup + ' --dntag ' + fdn + ' --metric wteslm --mkdir false'
# # # os.system(cmd)
# # subprocess.check_call(cmd)
# subprocess.check_call(["rum", "-w", "-q", "local", "query_tool", "--uptag", fup, "--dntag", fdn, "--metric", args.query, "--mkdir", "false"])


### analyze query ###
# load in results
frslt = '/xchip/cogs/hogstrom/analysis/scratch/Nov20/dose_analysis_tool.1353449771597/nov20/my_analysis.query_tool.2012112017162991/result_ESLM.COMBINED_n85x398050.gctx'
rslt = gct.GCT()
rslt.read(frslt)
rSigIds = rslt.get_rids()

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
#loop through each of the UNIQUE perts - graph ranks of query
pertSet = set(qPert)
for pert in pertSet:
	cmpd1 = pert
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
	qStr = qPertID[iPo[0]] #set pertID
	if len(qStr) >= 13:
		qStr = qStr[0:13] #shorten qPertID
	#run pymongo query
	CM = mutil.CMapMongo()
	#cmpdSigIds = CM.find({'pert_id':qStr},{'sig_id':True})
	cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
	if len(cmpdSigIds) < 2:
		print cmpd1 + ' has one or no instances in the cmap database'
		continue
	#loop through each dose
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
			np = 1000
			ntop = [x for x in i1 if x <= np]
			nPr = float(len(ntop))/(len(i1)) #percent of instances at the top of the list
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





# # # iteractve comands to analyze profiler file
# # #run from command line like this: python -m cProfile -o py_profiler.profile py_profiler_dose.py
# import pstats
# p = pstats.Stats('py_profiler2.profile')
# # p.sort_stats('calls','cumulative')
# p.sort_stats('cumulative','calls')
# p.sort_stats('time')
# p.print_stats(20)