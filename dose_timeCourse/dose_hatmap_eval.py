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

work_dir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/scratch'
fname = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
db = gct.GCT()
db.read(fname)


# instantiate a progress object
prog = progress.DeterminateProgressBar('Template Heatmaps')

# # read the data
# gcto = gct.GCT(args.res)
# gcto.read()
gcto = db
cids = gcto.get_cids()

# grab the cids from the file and mine dose information from them.  Find all of 
# the unique perts
# cids = gcto.get_gctx_cid(args.res)
pert_descs = gcto.get_column_meta('pert_desc')
doses = [float(x.split(':')[2]) for x in cids]
perts = [x.split(':')[1] for x in cids]
unique_perts = list(set(perts))

# grab the rid for use below
# rids = gcto.get_gctx_rid(args.res)
rids = gcto.get_rids()
LinCorrMtrx = numpy.zeros(shape=(978,len(unique_perts)))
num_perts = len(unique_perts)
for i,unique_pert in enumerate(unique_perts):
	prog.update('analyzing {0}'.format(unique_pert),i,num_perts)

	# grab the z-scores and doses for the current pert and sort the pairs
	# by dose. put the cid_inds in the same sorted order
	cid_inds = [j for j,x in enumerate(cids) if unique_pert in x]
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

		#save correlation data for all compounds to a matrix
		LinCorrMtrx[:,i] = linear_probe_corrs 


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


## plots ###
for j in range(LinCorrMtrx.shape[1]):
	col1 = LinCorrMtrx[:,j]
	col1.sort()
	plt.plot(col1,'.')
plt.show()


### pull some DMSO files ### 
qStr = 'DMSO'
CM = mutil.CMapMongo()
dmsoLst = CM.find({'pert_id':{'$regex':qStr},'cell_id':'PC3'},{'sig_id':True})
nDmsos = len(dmsoLst)
# #save dmso list
# dmsoOut = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/scratch/dmsoPC3_list.grp'
# with open(dmsoOut,'w') as f:
# 	for d in dmsoLst:
# 		f.write(d + '\n')
# load in affogato matrix data
affogato = gct.GCT()
fAff = '/xchip/cogs/data/build/affogato/affogato_r1_score_n398050x22268.gctx'
LMrid = db.get_rids()
affogato.read(fAff,cid=dmsoLst,rid=LMrid)
mtrx = affogato.matrix

# loop through and test the null distribution of expression-template correlations with different dose numbers
i=0
m = 20 #number of theoritical doses to test
LinCorrMtrx = numpy.zeros(shape=(978,m))
for n in range(1,m):
	#pull randomly from the DMSO
	# n = 4 #number of cids to pull
	iN = numpy.random.randint(0,nDmsos,n)
	nMtrx = mtrx[:,iN]

	# build prototype curves if there is more than one dose
	linear = numpy.linspace(1,10,n)
	log_gen = _log_gen(1)
	log_curve = [log_gen.next() for x in range(n)]
	log_gen = _log_gen(.5)
	half_log_curve = [log_gen.next() for x in range(n)]
	log_gen = _log_gen(.25)
	quarter_log_curve = [log_gen.next() for x in range(n)]
	curves = numpy.array([linear,log_curve,
						  half_log_curve,quarter_log_curve])
	num_probes = nMtrx.shape[0]
	cc = numpy.corrcoef(nMtrx,curves)

	# grab the correlation values for all the probes against prototype curves
	linear_probe_corrs = cc[0:num_probes,num_probes]
	log_probe_corrs = cc[0:num_probes,num_probes + 1]
	half_log_probe_corrs = cc[0:num_probes,num_probes + 2]
	quarter_log_probe_corrs = cc[0:num_probes,num_probes + 3]

	#save correlation data for all compounds to a matrix
	LinCorrMtrx[:,i] = linear_probe_corrs 
	i = i+1

# plt.interactive(True)
# for j in range(LinCorrMtrx.shape[1]):
# 	col1 = LinCorrMtrx[:,j]
# 	col1.sort()
# 	plt.plot(col1,'.')
# plt.show()

# plt.plot(LinCorrMtrx[:,1:4])
# plt.legend(['1','2','3'])
# plt.show()



n = 4 #number of cids to pull (the number of doses tested)
#calculate template curves
linear = numpy.linspace(1,10,n)
log_gen = _log_gen(1)
log_curve = [log_gen.next() for x in range(n)]
log_gen = _log_gen(.5)
half_log_curve = [log_gen.next() for x in range(n)]
log_gen = _log_gen(.25)
quarter_log_curve = [log_gen.next() for x in range(n)]
curves = numpy.array([linear,log_curve,
					  half_log_curve,quarter_log_curve])
num_probes = nMtrx.shape[0]

# loop through and test the null distribution many times with dose=4
perms = 1000 #number of theoritical doses to test
LinCorrMtrx = numpy.zeros(shape=(978,perms))
for p in range(0,perms):
	#pull randomly from the DMSO
	
	iN = numpy.random.randint(0,nDmsos,n)
	nMtrx = mtrx[:,iN]

	# build prototype curves if there is more than one dose
	cc = numpy.corrcoef(nMtrx,curves)

	# grab the correlation values for all the probes against prototype curves
	linear_probe_corrs = cc[0:num_probes,num_probes]
	log_probe_corrs = cc[0:num_probes,num_probes + 1]
	half_log_probe_corrs = cc[0:num_probes,num_probes + 2]
	quarter_log_probe_corrs = cc[0:num_probes,num_probes + 3]

	#save correlation data for all compounds to a matrix
	LinCorrMtrx[:,p] = linear_probe_corrs 

# for j in range(20,200):
for j in range(100):
	col1 = LinCorrMtrx[:,j]
	col1.sort()
	plt.plot(col1,)
plt.xlabel('probes ranked acording to corr with linear curve')
plt.ylabel('corr')
plt.show()

#corr threshold - at a give corr threshold, how many probes would you expect to pass by chance?
trVals = [.5, .7, .8, .9, .95]
threshMtx = numpy.zeros(shape=(LinCorrMtrx.shape[1],len(trVals))) 
for i,tr in enumerate(trVals):
	cThresh = tr
	trLst = []
	for j in range(LinCorrMtrx.shape[1]):
		col1 = LinCorrMtrx[:,j]
		col1.sort()
		colTr= col1[col1>cThresh]
		trLst.append(len(colTr))
	threshMtx[:,i] = numpy.array(trLst)
plt.hist(threshMtx[:,3])
plt.title('threshold correlation of .9 - DMSO simulation with 4 doses')
plt.xlabel('probes passing threshold')
plt.ylabel('freq from 1000 simulations')

