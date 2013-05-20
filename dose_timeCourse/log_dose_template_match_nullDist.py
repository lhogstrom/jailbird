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

cell = 'PC3'
tim = '6H'
cellLine = cell
timeP = tim
refControl = 'pc' #use pc vs vc controled data
# gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/by_pert_id_pert_dose/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
# load in zscore roast data
# gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/PRISM001_%s_%s_ZSPCQNORM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
# load in brewed by rna well
gctfile = glob.glob('/xchip/obelix/pod/brew_tmp/%s/PRISM001_%s_%s/by_rna_well/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
gctfile = gctfile[0]
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/ec50_test/%s_%s_%s' % (cell,timeP,refControl)
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

#permute the dmso lots of times
n_perm = 1000
n_curves = 4 #linear, log, half log, quarter log
template_names = ['linear', 'log', 'half_log', 'quarter_log']
dmso_perm_mtx = np.zeros((len(rids),n_perm,n_curves))
dmso_inds = [i for i,x in enumerate(perts) if 'DMSO' in x]
nRandNeeded = curves.shape[1]
for n in range(n_perm):
	iPermInd = dmso_inds
	np.random.shuffle(iPermInd) #shuffle indices of DMSOs
	iRandDMSO = iPermInd[0:nRandNeeded]
	dmso_data = gcto.matrix[:,iRandDMSO]
	dmso_cc = np.corrcoef(dmso_data,curves)
	linear_probe_corrs = dmso_cc[0:num_probes,num_probes]
	log_probe_corrs = dmso_cc[0:num_probes,num_probes + 1]
	half_log_probe_corrs = dmso_cc[0:num_probes,num_probes + 2]
	quarter_log_probe_corrs = dmso_cc[0:num_probes,num_probes + 3]
	dmso_perm_mtx[:,n,0] = linear_probe_corrs
	dmso_perm_mtx[:,n,1] = log_probe_corrs
	dmso_perm_mtx[:,n,2] = half_log_probe_corrs
	dmso_perm_mtx[:,n,3] = quarter_log_probe_corrs

# examineList = ['VX-680','MLN8054','GSK-1070916','BRD-K01737880','BRD-K29830875','AZD-1152HQPA']
examineList = ['BRD-K59369769-001-05-4','BRD-K83963101-001-01-0','BRD-K36740062-001-02-5','BRD-K01737880','BRD-K29830875','CMAP-AZD-1152HQPA']
num_perts = len(unique_perts)
doseLenList = []
observed_mtx = np.zeros((len(rids),len(examineList),n_curves))
for icmpd,cmpd in enumerate(examineList):
	ipert = unique_perts.index(cmpd)
	unique_pert = unique_perts[ipert]
	# for i,unique_pert in enumerate(unique_perts):
	# prog.update('analyzing {0}'.format(unique_pert),i,num_perts)
	# grab the z-scores and doses for the current pert and sort the pairs
	# by dose. put the cid_inds in the same sorted order
	cid_inds = [i for i,x in enumerate(perts) if unique_pert in x]
	doseLength = len(cid_inds)
	doseLenList.append(doseLength)
	pert_desc = pert_descs[cid_inds[0]] #set pert desc to the first dose
	pert_doses = [float(doses[x]) for x in cid_inds]
	tmp_tup = zip(pert_doses,cid_inds)
	tmp_tup.sort() #is this sorting strings? - not correct order of doses
	pert_doses,cid_inds = zip(*tmp_tup)
	# if len(pert_doses) > 1:
	# build prototype curves if there is more than one dose
	linear = np.linspace(1,10,len(pert_doses))
	log_gen = _log_gen(1)
	log_curve = [log_gen.next() for x in range(len(pert_doses))]
	log_gen = _log_gen(.5)
	half_log_curve = [log_gen.next() for x in range(len(pert_doses))]
	log_gen = _log_gen(.25)
	quarter_log_curve = [log_gen.next() for x in range(len(pert_doses))]
	curves = np.array([linear,log_curve,
						  half_log_curve,quarter_log_curve])
	# curves = np.array([linear])
	# correlate all of the probes in the data to the prototype curves
	pert_data = gcto.matrix[:,cid_inds]
	num_probes = pert_data.shape[0]
	cc = np.corrcoef(pert_data,curves) #does spearman or pearson make more sense here?
	linear_probe_corrs = cc[0:num_probes,num_probes]
	log_probe_corrs = cc[0:num_probes,num_probes + 1]
	half_log_probe_corrs = cc[0:num_probes,num_probes + 2]
	quarter_log_probe_corrs = cc[0:num_probes,num_probes + 3]
	observed_mtx[:,icmpd,0] = linear_probe_corrs
	observed_mtx[:,icmpd,1] = log_probe_corrs
	observed_mtx[:,icmpd,2] = half_log_probe_corrs
	observed_mtx[:,icmpd,3] = quarter_log_probe_corrs
	### calculate p and fdr 
	p_vec = np.zeros_like(linear_probe_corrs)
	for i,rho in enumerate(linear_probe_corrs):
		if rho > 0:
			pass_obs_mtrx = np.greater(dmso_perm_mtx,rho)
			p = (np.sum(pass_obs_mtrx)) / float(pass_obs_mtrx.size)*2
		if rho < 0:
			pass_obs_mtrx = np.less(dmso_perm_mtx,rho)
			p = (np.sum(pass_obs_mtrx)) / float(pass_obs_mtrx.size)*2
		p_vec[i] = p
	q = .1 #FDR threshold	
	pID, pN = FDR(p_vec,q)
	#	# The two-sided p-value of the test is calculated as the proportion of sampled 
	#   # permutations where the absolute difference was greater than or equal to ABS(T(obs))
		#single -sided test - p = number of rands that pass rho / total # rands
	#compare empirical p-values to corr p-values


	#loop through each of the templates - plot real vs null dist
	for n_templ in range(n_curves):
		templates = template_names[n_templ]
		randFlat = dmso_perm_mtx[:,:,n_templ].flatten()
		ObsFlat = observed_mtx[:,icmpd,n_templ].flatten()
		maxRand = max(randFlat)
		maxObs = max(ObsFlat)
		minRand = min(randFlat)
		minObs = min(ObsFlat)
		range1=[min([minObs,minRand]), max([maxObs,maxRand])]
		fig = plt.figure()
		ax = fig.add_subplot(2,1,1)
		weightObs = np.zeros_like(ObsFlat) + 1. / len(ObsFlat)
		plt.hist(ObsFlat,bins=30,range=range1,weights=weightObs)
		plt.xlabel('observed es scores')
		# n, bins, patches = plt.hist(colCMAPscore,bins=30,range=range1,weights=weightObs)
		# plt.bar(bins,n)
		ax = fig.add_subplot(2,1,2)
		weightRand = np.zeros_like(randFlat) + 1. / len(randFlat)
		plt.hist(randFlat,bins=30,range=range1,weights=weightRand)
		plt.xlabel('random dmso')
		work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/log_template_fit'
		plt.savefig(os.path.join(work_dir,cmpd + '_' + templates + '_templates_cor.png'))
		plt.close
	# plt.show()








### scratch work #### 
#corr rank for first compound
isort = numpy.argsort(linear_probe_corrs)[::-1]
n_show = 15
itop = isort[:n_show]
col1 = pert_data[itop[1],:]
topData = pert_data[itop,:]
n_templ = 0
scipy.stats.pearsonr(list(col1),list(curves[n_templ]))
templates = template_names[n_templ]
plt.plot(pert_doses,col1,'o')
plt.show()
topRids = [rids[i] for i in itop]
heatStruc = gct.GCT()
heatStruc.build(topData,topRids,pert_doses,{},{})
heatStructF = os.path.join(work_dir,cmpd + '_' + templates + '_heatmap_data.gctx')
heatStruc.write(heatStructF)
outHeatmap = os.path.join(work_dir,cmpd + '_' + templates + '_heatmap')
heatmap_cmd = 'heatmap -o ' + outHeatmap + ' --linkage 0 --gct ' + heatStructF + ' --column-text id --row-text id --format png'
#/xchip/cogs/bin/heatmap --gct /xchip/cogs/hogstrom/analysis/scratch/prism/ec50_test/PC3_6H_pc/BRD-K59369769-001-05-4_linear_heatmap_data.gctx -o /xchip/cogs/hogstrom/analysis/scratch/prism/ec50_test/PC3_6H_pc/BRD-K59369769-001-05-4_linear_heatmap --column-text id --row-text id --format png
#make plots of raw data with templates
for icmpd,cmpd in enumerate(examineList):
	for n_templ in range(n_curves):
		corrValues = observed_mtx[:,icmpd,n_templ]
		isort = numpy.argsort(corrValues)[::-1]
		n_show = 5
		itop = isort[:n_show]
		#get raw data
		ipert = unique_perts.index(cmpd)
		unique_pert = unique_perts[ipert]
		cid_inds = [i for i,x in enumerate(perts) if unique_pert in x]
		pert_data = gcto.matrix[:,cid_inds]
		pert_doses = [float(doses[x]) for x in cid_inds]
		for igene in itop:
			obsZs = pert_data[igene,:]
			# cc = np.corrcoef(obsZs,pert_doses)
			# cc = scipy.stats.pearsonr(obsZs,pert_doses)
			cc = scipy.stats.pearsonr(obsZs,list(curves[n_templ,:]))

### write p_vec for matlab fdr comparsion ##
pvecFile = os.path.join(work_dir,'pval_vec.txt')
with open(pvecFile,'a') as f:
	for p in p_vec:
		f.write(str(p) + '\n')
	f.write('\n')

### qq plot
logpO = -np.log10(np.sort(p_vec))
nullP = np.array(range(1,len(p_vec)+1))/float(len(p_vec))
logpE = -np.log10(nullP)


# # http://gettinggeneticsdone.blogspot.com/2010/07/qq-plots-of-p-values-in-r-using-base.html
# # Define the function
# ggd.qqplot = function(pvector, main=NULL, ...) {
#     o = -log10(sort(pvector,decreasing=F))
#     e = -log10( 1:length(o)/length(o) )
#     plot(e,o,pch=19,cex=1, main=main, ...,
#         xlab=expression(Expected~~-log[10](italic(p))),
#         ylab=expression(Observed~~-log[10](italic(p))),
#         xlim=c(0,max(e)), ylim=c(0,max(o)))
#     lines(e,e,col="red")
# }
 
# #Generate some fake data that deviates from the null
# set.seed(42)
# pvalues=runif(10000)
# pvalues[sample(10000,10)]=pvalues[sample(10000,10)]/5000
 
# # pvalues is a numeric vector
# pvalues[1:10]
 
# # Using the ggd.qqplot() function
# ggd.qqplot(pvalues)
 
# # Add a title
# ggd.qqplot(pvalues, "QQ-plot of p-values using ggd.qqplot")


