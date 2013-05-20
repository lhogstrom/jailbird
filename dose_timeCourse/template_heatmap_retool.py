#! /usr/bin/env python
'''
run dose-template matching - test for null distribution
'''
import os
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import cmap.io.gct as gct
import cmap.util.tool_ops as tool_ops
import cmap.util.progress as progress
import cmap.analytics.fitting as fitting
import cmap.analytics.cluster as cluster
import cmap.util.queue as queue
import cmap.io.plategrp as grp
import cmap.plot.colors as colors
import random
import scipy
### template matching - null distirbution
#plate data - colapsed by rna well
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

prog = progress.DeterminateProgressBar('Template Heatmaps')

cellLst = ['PC3', 'A375', 'MCF7']
timeLst = ['6H', '24H']
for cell in cellLst:
	for tim in timeLst:
		# cell = 'PC3'
		# tim = '6H'
		cellLine = cell
		timeP = tim
		#load data
		refControl = 'pc' #use pc vs vc controled data
		gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/by_pert_id_pert_dose/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		# load in zscore roast data
		# gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/PRISM001_%s_%s_ZSPCQNORM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		# load in brewed by rna well
		# gctfile = glob.glob('/xchip/obelix/pod/brew_tmp/%s/PRISM001_%s_%s/by_rna_well/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		# load in brewed by rna well - INFERED
		# gctfile = glob.glob('/xchip/obelix/pod/brew_tmp/%s/PRISM001_%s_%s/by_rna_well/PRISM001_%s_%s_COMPZ.MODZ_SCORE_n*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))		
		gctfile = gctfile[0]
		work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/heatmap_by_pert_id_pert_dose/%s_%s_%s' % (cell,timeP,refControl)
		if not os.path.exists(work_dir):
			os.mkdir(work_dir)
		gcto = gct.GCT() #make a gct object
		gcto.read(gctfile)

		cids = gcto.get_cids()
		pert_descs = gcto.get_column_meta('pert_desc')
		doses = gcto.get_column_meta('pert_dose')
		perts = gcto.get_column_meta('pert_id')
		# doses = [float(x.split(':')[2]) for x in cids]
		# perts = [x.split(':')[1] for x in cids]
		unique_perts = list(set(perts))

		# grab the rid for use below
		rids = gcto.get_rids()

		### changed to accomodate different plate structure
		examineList = ['VX-680','MLN8054','GSK-1070916','BRD-K01737880','BRD-K29830875','AZD-1152HQPA']
		# examineList = ['BRD-K59369769-001-05-4','BRD-K83963101-001-01-0','BRD-K36740062-001-02-5','BRD-K01737880','BRD-K29830875','CMAP-AZD-1152HQPA']
		# examineList = ['BRD-K36740062-001-02-5']
		num_perts = len(examineList)
		# template_names = ['linear', 'log', 'half_log', 'quarter_log']
		# template_names = ['linear', 'log10']
		# doseLenList = []
		# observed_mtx = np.zeros((len(rids),len(examineList),len(template_names)))
		for icmpd,unique_pert in enumerate(examineList):
			prog.update('analyzing {0}'.format(unique_pert),icmpd,num_perts)
			
			#### changes to the code are bellow this point ####
			
			# grab the z-scores and doses for the current pert and sort the pairs
			# by dose. put the cid_inds in the same sorted order
			# cid_inds = [i for i,x in enumerate(cids) if unique_pert in x]
			cid_inds = [i for i,x in enumerate(perts) if unique_pert in x]
			# pert_desc = pert_descs[cid_inds[0]] #set pert desc to the first dose
			pert_doses = [float(doses[x]) for x in cid_inds]
			tmp_tup = zip(pert_doses,cid_inds)
			tmp_tup.sort()
			pert_doses,cid_inds = zip(*tmp_tup)
			pert_data = gcto.matrix[:,cid_inds]
			if len(pert_doses) > 1:
				# build prototype curves if there is more than one dose
				# template_names = ['linear', 'log', 'half_log', 'quarter_log']
				# log_steps = [0, 1, .5, .25] #set to 0 if linear
				# for istep,step in enumerate(log_steps):
				template_names = ['linear', 'log10', 'log2']
				# template_names = ['linear']
				for istep,step in enumerate(template_names):
					template1 = step
					if step == 'linear':
						template_curve = np.array(pert_doses)
					elif step == 'log10':
						# log_gen = _log_gen(step)
						# template_curve = [log_gen.next() for x in range(len(pert_doses))]
						template_curve = np.log10(pert_doses)
					elif step == 'log2':
						template_curve = np.log2(pert_doses)
					else:
						print 'template name error'
					cc_list = [scipy.stats.pearsonr(pert_data[x,:],template_curve) for x in range(len(rids))]
					rho_vec = [cc_list[x][0] for x in range(len(rids))]
					rho_vec = np.array(rho_vec)
					p_vec = [cc_list[x][1] for x in range(len(rids))]
					p_vec = np.array(p_vec)
					q = .1 #FDR threshold
					pID, pN = FDR(p_vec,q) #find FDR threshold
					if type(pID) == list:
						print unique_pert + 'matching to ' + template1 + ' template - perterbation does not have any significant genes that pass the FDR threshold'
						continue
					else:
						pass_fdr = np.less_equal(p_vec,pID) 
						ipass_fdr = np.array(range(len(rids)))[pass_fdr] #get indices which pass fdr
						iRhoSort = np.argsort(rho_vec[ipass_fdr])[::-1]
						iRhoSorted_passFDR = ipass_fdr[iRhoSort] #these are indices which pass FDR and are sorted by correlation
						data_pass_fdr = pert_data[iRhoSorted_passFDR,:]
						ordered_rids = [rids[i] for i in iRhoSorted_passFDR]
						heatStruc = gct.GCT()
						heatStruc.build(data_pass_fdr,ordered_rids,pert_doses,{},{})
						heatStructF = os.path.join(work_dir,unique_pert + '_' + template1 + '_heatmap_data.gctx')
						heatStruc.write(heatStructF)
						outHeatmap = os.path.join(work_dir,unique_pert + '_' + template1 + '_heatmap')
						heatmap_cmd = 'heatmap -o ' + outHeatmap + ' --linkage 0 --gct ' + heatStructF + ' --column-text id --row-text id --format png'
						os.system(heatmap_cmd)
						# # write the p values and correlations out to file
						with open(os.path.join(work_dir,unique_pert +'_'+ template1 + '_template_match_summary.txt'),'w') as f:
							f.write('\t'.join(['probeset','template corr', 'tempalte p', 'fdir q = ' + str(q)]) + '\n')
							for j in iRhoSorted_passFDR:
								f.write('\t'.join([rids[j],str(rho_vec[j]), str(p_vec[j])]) + '\n')

# # loop through and read in stats about dose responsive genes
# cellLst = ['PC3', 'A375', 'MCF7']
# timeLst = ['6H', '24H']
# examineList = ['BRD-K59369769-001-05-4','BRD-K83963101-001-01-0','BRD-K36740062-001-02-5','BRD-K01737880','BRD-K29830875','CMAP-AZD-1152HQPA']
# upProbeDict = {}
# for icmpd,unique_pert in enumerate(examineList):
# 	upProbeDict[unique_pert] = {}
# for cell in cellLst:
# 	for tim in timeLst:
# 		# cell = 'PC3'
# 		# tim = '6H'
# 		cellLine = cell
# 		timeP = tim
# 		refControl = 'pc' 
# 		work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/ec50_test/%s_%s_%s' % (cell,timeP,refControl)
# 		template_names = ['linear', 'log', 'half_log', 'quarter_log']
# 		for icmpd,unique_pert in enumerate(examineList):
# 			template1 = 'linear'
# 			sumFile = os.path.join(work_dir,unique_pert +'_'+ template1 + '_template_match_summary.txt')
# 			probes = []
# 			corrs = []
# 			pvals = []
# 			if os.path.isfile(sumFile): #files were created if genes pass fdr
# 				with open(sumFile) as f:
# 					for i,row in enumerate(f):
# 						if i == 0:
# 							continue
# 						else:
# 							probes.append(row.split('\t')[0]) 
# 							corrs.append(float(row.split('\t')[1]))
# 							pvals.append(float(row.split('\t')[2][:-1]))
# 			probesUp = [x for i,x in enumerate(probes) if corrs[i] > 0] #probes positivly correlating with template
# 			probesDn = [x for i,x in enumerate(probes) if corrs[i] < 0]
# 			# nestedDict = {}
# 			# nestedDict[cell+'_'+tim] = probesUp
# 			upProbeDict[unique_pert][cell+'_'+tim] = probesUp

# for cmpd in upProbeDict:
# 	for context in upProbeDict[cmpd]:
# 		probeList = upProbeDict[cmpd][context]


# keys1 = upProbeDict[cmpd].keys()
# set1 = set(upProbeDict[cmpd][keys1[0]])
# set2 = set(upProbeDict[cmpd][keys1[1]])
# set3 = set(upProbeDict[cmpd][keys1[3]])





		# 	rho_pass_fdr = rho_vec[pass_fdr]
		# 	data_pass_fdr = pert_data[pass_fdr,:]
		# 	iRhoSort = np.argsort(rho_vec)[::-1]
		# 	iRhoSorted_passFDR = iRhoSort[pass_fdr]
		# 	rho_vec[iRhoSorted_passFDR]

		# linear_probe_corrs_sort_ind = numpy.argsort(linear_probe_corrs_p)[::-1]
		# top = pert_data[linear_probe_corrs_sort_ind[0:50],:]
		# bot = pert_data[linear_probe_corrs_sort_ind[-50:],:]


		# curves = numpy.array([linear,log_curve,
		# 					  half_log_curve,quarter_log_curve])

		# # correlate all of the probes in the data to the prototype curves
		# pert_data = gcto.matrix[:,cid_inds]
		# num_probes = pert_data.shape[0]
		# cc = numpy.corrcoef(pert_data,curves)
		# scipy.stats.pearsonr(pert_data,curves)
		# # grab the correlation values for all the probes against prototype curves
		# linear_probe_corrs = cc[0:num_probes,num_probes]
		# log_probe_corrs = cc[0:num_probes,num_probes + 1]
		# half_log_probe_corrs = cc[0:num_probes,num_probes + 2]
		# quarter_log_probe_corrs = cc[0:num_probes,num_probes + 3]

		# # compute the random correlation profile for this pert
		# num_probes = gcto.matrix.shape[0]
		# probe_inds = range(num_probes)
		# linear_perm_cc = []
		# log_perm_cc = []
		# half_log_perm_cc = []
		# quarter_log_perm_cc = []
		# for i in range(1000):
		# 	perm_curve_inds = [random.sample(probe_inds,1)[0] for x in range(len(pert_doses))]
		# 	perm_curve = [pert_data[perm_curve_inds[x],x] for x in range(len(pert_doses))]
		# 	perm_covar = numpy.corrcoef(perm_curve,curves)
		# 	linear_perm_cc.append(perm_covar[0][1])
		# 	log_perm_cc.append(perm_covar[0][2])
		# 	half_log_perm_cc.append(perm_covar[0][3])
		# 	quarter_log_perm_cc.append(perm_covar[0][4])

		# # compute the nominal p values for all correlation values
		# linear_probe_corrs_p = numpy.array([stats.percentileofscore(linear_perm_cc,x) 
		# 						for x in linear_probe_corrs])
		# log_probe_corrs_p = numpy.array([stats.percentileofscore(log_perm_cc,x) 
		# 						for x in log_probe_corrs])
		# half_log_probe_corrs_p = numpy.array([stats.percentileofscore(half_log_perm_cc,x) 
		# 						for x in half_log_probe_corrs])
		# quarter_log_probe_corrs_p = numpy.array([stats.percentileofscore(quarter_log_perm_cc,x) 
		# 						for x in quarter_log_probe_corrs])

		# # write the p values and correlations out to file
		# with open(os.path.join(work_dir,pert_desc + '_template_match_summary.txt'),'w') as f:
		# 	f.write('\t'.join(['probeset','linear corr', 'linear p','log corr', 'log p',
		# 		'half-log corr', 'half-log p','quarter-log corr', 'quarter-log p']) + '\n')
		# 	for j in range(len(linear_probe_corrs)):
		# 		f.write('\t'.join([rids[j],str(linear_probe_corrs[j]), str(linear_probe_corrs_p[j])
		# 			,str(log_probe_corrs[j]), str(log_probe_corrs_p[j])
		# 			,str(half_log_probe_corrs[j]), str(half_log_probe_corrs_p[j])
		# 			,str(quarter_log_probe_corrs[j]), str(quarter_log_probe_corrs_p[j])]) + '\n')


		# # build the linear heatmap
		# linear_probe_corrs_sort_ind = numpy.argsort(linear_probe_corrs_p)[::-1]
		# top = pert_data[linear_probe_corrs_sort_ind[0:50],:]
		# bot = pert_data[linear_probe_corrs_sort_ind[-50:],:]
		# combined = numpy.vstack([top,bot])
		# combined_row_normalized =  combined + numpy.abs(numpy.array([numpy.min(combined,1)]).T)
		# row_sums = combined_row_normalized.sum(axis=1)
		# combined_row_normalized =  combined_row_normalized / row_sums[:,numpy.newaxis]
		# plt.imshow(combined_row_normalized,interpolation='nearest',cmap='RdBu')
		# plt.axis('off')
		# plt.savefig(os.path.join(work_dir,pert_desc + '_linear_heatmap.png'))

		# # build the log heatmap
		# log_probe_corrs_sort_ind = numpy.argsort(log_probe_corrs_p)[::-1]
		# top = pert_data[log_probe_corrs_sort_ind[0:50],:]
		# bot = pert_data[log_probe_corrs_sort_ind[-50:],:]
		# combined = numpy.vstack([top,bot])
		# combined_row_normalized =  combined + numpy.abs(numpy.array([numpy.min(combined,1)]).T)
		# row_sums = combined_row_normalized.sum(axis=1)
		# combined_row_normalized =  combined_row_normalized / row_sums[:,numpy.newaxis]
		# plt.imshow(combined_row_normalized,interpolation='nearest',cmap='RdBu')
		# plt.axis('off')
		# plt.savefig(os.path.join(work_dir,pert_desc + '_log_heatmap.png'))

		# # build the half log heatmap
		# half_log_probe_corrs_sort_ind = numpy.argsort(half_log_probe_corrs_p)[::-1]
		# top = pert_data[half_log_probe_corrs_sort_ind[0:50],:]
		# bot = pert_data[half_log_probe_corrs_sort_ind[-50:],:]
		# combined = numpy.vstack([top,bot])
		# combined_row_normalized =  combined + numpy.abs(numpy.array([numpy.min(combined,1)]).T)
		# row_sums = combined_row_normalized.sum(axis=1)
		# combined_row_normalized =  combined_row_normalized / row_sums[:,numpy.newaxis]
		# plt.imshow(combined_row_normalized,interpolation='nearest',cmap='RdBu')
		# plt.axis('off')
		# plt.savefig(os.path.join(work_dir,pert_desc + '_half_log_heatmap.png'))

		# # build the quarter log heatmap
		# quarter_log_probe_corrs_sort_ind = numpy.argsort(quarter_log_probe_corrs_p)[::-1]
		# top = pert_data[quarter_log_probe_corrs_sort_ind[0:50],:]
		# bot = pert_data[quarter_log_probe_corrs_sort_ind[-50:],:]
		# combined = numpy.vstack([top,bot])
		# combined_row_normalized =  combined + numpy.abs(numpy.array([numpy.min(combined,1)]).T)
		# row_sums = combined_row_normalized.sum(axis=1)
		# combined_row_normalized =  combined_row_normalized / row_sums[:,numpy.newaxis]
		# plt.imshow(combined_row_normalized,interpolation='nearest',cmap='RdBu')
		# plt.axis('off')
		# plt.savefig(os.path.join(work_dir,pert_desc + '_quarter_log_heatmap.png'))

		# # clear that progress bar
		# prog.clear()

