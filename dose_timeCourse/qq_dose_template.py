#! /usr/bin/env python
'''
make a qq plot where the dot size is controled by the dose 
correlation with a linear template
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
import cmap.analytics.dose as doseClass
import cmap.analytics.cluster as cluster
import cmap.util.queue as queue
import cmap.io.plategrp as grp
import cmap.plot.colors as colors
import random
import scipy
import scipy.stats as stats

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

def FDR(p,q):
	'''
	input:
	p   - vector of p-values (numpy array)
	q   - False Discovery Rate level
	output:
	pID - p-value threshold based on independence or positive dependence
	pN  - Nonparametric p-value threshold
	*if non of the p-values pass FDR, an empty list is returned
	'''
	isort = np.argsort(p)
	p = p[isort]
	V = len(p)
	I = np.array(range(1,V+1))
	cVID = 1
	cVN = sum(np.divide(float(1),I))
	#threshold based on independence or positive dependence
	ID_thresh_vec = (I*q)/V/cVID
	ID_pass_vec = np.less_equal(p,ID_thresh_vec)
	if any(ID_pass_vec):
		pID = max(p[ID_pass_vec])
	else:
		pID = []
	# Nonparametric threshold
	N_thresh_vec = (I*q)/V/cVN
	N_pass_vec = np.less_equal(p,N_thresh_vec)
	if any(N_pass_vec):
		pN = max(p[N_pass_vec])
	else:
		pN = []
	return pID, pN

prog = progress.DeterminateProgressBar('Template Heatmaps')


cell = 'PC3'
tim = '6H'
cellLine = cell
timeP = tim
#load data
refControl = 'pc' #use pc vs vc controled data
# gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/by_pert_id_pert_dose/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
# load in zscore roast data
# gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/PRISM001_%s_%s_ZSPCQNORM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
# load in brewed by rna well
gctfile = glob.glob('/xchip/obelix/pod/brew_tmp/%s/PRISM001_%s_%s/by_rna_well/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
# load in brewed by rna well - INFERED
# gctfile = glob.glob('/xchip/obelix/pod/brew_tmp/%s/PRISM001_%s_%s/by_rna_well/PRISM001_%s_%s_COMPZ.MODZ_SCORE_n*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))		
gctfile = gctfile[0]
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/ec50_test_inf/%s_%s_%s' % (cell,timeP,refControl)
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
# examineList = ['VX-680','MLN8054','GSK-1070916','BRD-K01737880','BRD-K29830875','AZD-1152HQPA']
examineList = ['BRD-K59369769-001-05-4','BRD-K83963101-001-01-0','BRD-K36740062-001-02-5','BRD-K01737880','BRD-K29830875','CMAP-AZD-1152HQPA']
# examineList = ['BRD-K36740062-001-02-5']
num_perts = len(examineList)
template_names = ['linear', 'log', 'half_log', 'quarter_log']
doseLenList = []
observed_mtx = np.zeros((len(rids),len(examineList),len(template_names)))

icmpd = 0
unique_pert = examineList[icmpd]
cid_inds = [i for i,x in enumerate(perts) if unique_pert in x]
pert_desc = pert_descs[cid_inds[0]] #set pert desc to the first dose
pert_doses = [float(doses[x]) for x in cid_inds]
tmp_tup = zip(pert_doses,cid_inds)
tmp_tup.sort()
pert_doses,cid_inds = zip(*tmp_tup)
pert_data = gcto.matrix[:,cid_inds]
# if len(pert_doses) > 1:
# 	# build prototype curves if there is more than one dose
template_names = ['linear', 'log', 'half_log', 'quarter_log']
log_steps = [0, 1, .5, .25] #set to 0 if linear
# for istep,step in enumerate(log_steps):
istep = 0
step = log_steps[istep]
template1 = template_names[istep]
if step == 0:
	template_curve = np.linspace(1,10,len(pert_doses))
else:
	log_gen = _log_gen(step)
	template_curve = [log_gen.next() for x in range(len(pert_doses))]
cc_list = [scipy.stats.pearsonr(pert_data[x,:],template_curve) for x in range(len(rids))]
rho_vec = [cc_list[x][0] for x in range(len(rids))]
rho_vec = np.array(rho_vec)
p_vec = [cc_list[x][1] for x in range(len(rids))]
p_vec = np.array(p_vec)
q = .1 #FDR threshold
pID, pN = FDR(p_vec,q) #find FDR threshold
if type(pID) == list:
	print unique_pert + 'matching to ' + template1 + ' template - perterbation does not have any significant genes that pass the FDR threshold'
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
	# os.system(heatmap_cmd)
	# # # write the p values and correlations out to file
	# with open(os.path.join(work_dir,unique_pert +'_'+ template1 + '_template_match_summary.txt'),'w') as f:
	# 	f.write('\t'.join(['probeset','template corr', 'tempalte p', 'fdir q = ' + str(q)]) + '\n')
	# 	for j in iRhoSorted_passFDR:
	# 		f.write('\t'.join([rids[j],str(rho_vec[j]), str(p_vec[j])]) + '\n')

zVec = pert_data[:,-1] #set zVec to be the highest dose
iSort = np.argsort(zVec)
zVec = zVec[iSort]
mu, sigma = 0, 1
nullDist = np.random.normal(mu, sigma, len(rids))
nullDist.sort()

y_label = 'observed z-score'
x_label = 'expected z-score'
legend_string1 = 'positive dose correlation'
legend_string2 = 'negative dose correlation'
graph_title = unique_pert # + ' dose = ' + dose1
# outf = os.path.join(work_dir,fname)
#set all negative corr values to 0
rho_vec_zSort = rho_vec[iSort] # sort rho values by index of z-score sort
posRho = np.zeros_like(rho_vec_zSort)
iposRho = np.greater_equal(rho_vec_zSort,0)
posRho[iposRho] = rho_vec_zSort[iposRho]
#set all positive corr values to 0
negRho = np.zeros_like(rho_vec_zSort)
inegRho = np.less_equal(rho_vec_zSort,0)
negRho[inegRho] = rho_vec_zSort[inegRho]
#set markers size to be based on correlation with template
marker_size1 = np.multiply(posRho,100)
marker_size2 = np.multiply(np.absolute(negRho),100)

import cmap.plot.qq as qq
qq.qq_plot(zVec,marker_size1,marker_size2,nullDist,graph_title,y_label, x_label, legend_string1,legend_string2)



### scratch for qq stuff ### 
cids = gcto.get_cids()
pert_descs = gcto.get_column_meta('pert_desc')
doses = gcto.get_column_meta('pert_dose')
perts = gcto.get_column_meta('pert_id')
unique_perts = list(set(perts))
rids = gcto.get_rids()

dp = doseClass.DosePlate()
dp.add_from_gct(gcto.src)
dp.examine_doses_tested()
examineList = dp.perts_at_dose
num_perts = len(examineList)
for icmpd,unique_pert in enumerate(examineList):
	icmpd = 0
	unique_pert = examineList[icmpd]
	prog.update('qq plot {0}'.format(unique_pert),icmpd,num_perts)
	cid_inds = [i for i,x in enumerate(perts) if unique_pert in x]
	pert_desc = pert_descs[cid_inds[0]] #set pert desc to the first dose
	pert_doses = [float(doses[x]) for x in cid_inds]
	tmp_tup = zip(pert_doses,cid_inds)
	tmp_tup.sort()
	pert_doses,cid_inds = zip(*tmp_tup)
	pert_data = gcto.matrix[:,cid_inds]
	if len(pert_doses) > 1:
		template_names = ['linear', 'log10', 'log2']
		for istep,step in enumerate(template_names):
			template1 = step
			if step == 'linear':
				template_curve = np.array(pert_doses)
			elif step == 'log10':
				template_curve = np.log10(pert_doses)
			elif step == 'log2':
				template_curve = np.log2(pert_doses)
			else:
				print 'template name error'
			# calcualte stats on observation of interest 
			cc_list = [stats.pearsonr(pert_data[x,:],template_curve) for x in range(len(rids))]
			rho_vec = [cc_list[x][0] for x in range(len(rids))]
			rho_vec = np.array(rho_vec)
			p_vec = [cc_list[x][1] for x in range(len(rids))]
			p_vec = np.array(p_vec)
			# ###run permutations to creat null distribution of corr values
			nPerm = 1000
			permRhoMtrx = np.zeros((len(rids),nPerm))
			for perm in range(nPerm):
				iRandObs = range(pert_data.shape[1])
				np.random.shuffle(iRandObs)
				corrs = np.corrcoef(template_curve,pert_data[:,iRandObs])
				permRhoMtrx[:,perm] = corrs[0,1:]
				#test to see if two calculations methods are the same to a given precision
				# cc_list = [stats.pearsonr(pert_data[x,iRandObs],template_curve) for x in range(len(rids))] #this takes too long
				# rho_vec = [cc_list[x][0] for x in range(len(rids))]					
				# np.allclose(perm_list, np.array(rho_vec),rtol=1e-06)
			# ### calculate p-value based on null distribution
			# #compare observed gene to null gene 
			# grtrMtrx = np.greater(np.abs(permRhoMtrx.T),np.abs(rho_vec))
			# null_pVec1 = (1 + np.sum(grtrMtrx, axis=0)) / float(nPerm)
			# #compare observed gene to all null genes
			# null_pVec = np.zeros_like(rho_vec)
			# for igene in range(len(rho_vec)):
			# 	rho = rho_vec[igene]
			# 	p = np.sum(np.abs(permRhoMtrx.flatten()) > np.abs(rho)) /float(len(permRhoMtrx.flatten()))
			# 	null_pVec[igene] = p
			q = .1 #FDR threshold
			pID, pN = FDR(p_vec,q) #find FDR threshold
			#plot histogram of template correlations with null distribution
			n, bins, patches = plt.hist(permRhoMtrx.flatten(),60,normed=1)
			plt.close()
			plt.hist(rho_vec,40,normed=1)
			r1 = plt.plot(bins[1:],n,'r-',linewidth=3)
			plt.legend([r1], ['null distribution'], numpoints=1, loc=1)
			plt.title(unique_pert + ' ' + step)
			plt.xlabel('correltion with template')
			plt.ylabel('relative frequency')
			plt.savefig(os.path.join(work_dir,unique_pert + '_' + step + '_corr_dist.png'))
			plt.close()
			if step == 'linear':
				zVec = pert_data[:,-1] #set zVec to be the highest dose
				iSort = np.argsort(zVec)
				zVec = zVec[iSort]
				#make null distribution
				mu, sigma = 0, 1
				nullDist = np.random.normal(mu, sigma, len(rids))
				nullDist.sort()
				#set plotting inputs
				y_label = 'observed z-score'
				x_label = 'expected z-score'
				legend_string1 = 'positive dose correlation'
				legend_string2 = 'negative dose correlation'
				graph_title = unique_pert # + ' dose = ' + dose1
				outf = os.path.join(work_dir,unique_pert + 'qq_linear_dose.png')
				#set all negative corr values to 0
				rho_vec_zSort = rho_vec[iSort] # sort rho values by index of z-score sort
				posRho = np.zeros_like(rho_vec_zSort)
				iposRho = np.greater_equal(rho_vec_zSort,0)
				posRho[iposRho] = rho_vec_zSort[iposRho]
				#set all positive corr values to 0
				negRho = np.zeros_like(rho_vec_zSort)
				inegRho = np.less_equal(rho_vec_zSort,0)
				negRho[inegRho] = rho_vec_zSort[inegRho]
				#set markers size to be based on correlation with template
				marker_size1 = np.multiply(posRho,100)
				marker_size2 = np.multiply(np.absolute(negRho),100)
				qq.qq_plot(zVec,marker_size1,marker_size2,nullDist,graph_title,
					y_label, x_label, legend_string1,legend_string2,outf,showfig=False)
				plt.close()
				
				#qq of null distribution
				rho_vec_zSort = permRhoMtrx[iSort,0] # pick specific permutation
				posRho = np.zeros_like(rho_vec_zSort)
				iposRho = np.greater_equal(rho_vec_zSort,0)
				posRho[iposRho] = rho_vec_zSort[iposRho]
				#set all positive corr values to 0
				negRho = np.zeros_like(rho_vec_zSort)
				inegRho = np.less_equal(rho_vec_zSort,0)
				negRho[inegRho] = rho_vec_zSort[inegRho]
				#set markers size to be based on correlation with template
				marker_size1 = np.multiply(posRho,100)
				marker_size2 = np.multiply(np.absolute(negRho),100)
				qq.qq_plot(zVec,marker_size1,marker_size2,nullDist,graph_title,
					y_label, x_label, legend_string1,legend_string2)
				plt.close()


