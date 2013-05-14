#!/bin/py 
import glob
import os
import cmap.analytics.dose as doseClass
import cmap.io.gct as gct
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

tim = '6H'
cell = 'A375'
cellLine = cell
timeP = tim
refControl = 'vc' #use pc vs vc controled data
#load in data by rna well:
gctfile = glob.glob('/xchip/cogs/projects/PRISM/data_by_rna_well/%s/PRISM001_%s_%s/PRISM001_%s_%s_ZSVCQNORM_n*x978.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
#load in data by pert_id_pert_desc:
# gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/by_pert_id_pert_dose/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
gctfile = gctfile[0]
# gctfile = '/xchip/obelix/pod/brew/vc/PRISM001_A375_24H/by_pert_id_pert_dose/PRISM001_A375_24H_COMPZ.MODZ_SCORE_LM_n58x978.gctx'
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/%s_%s_%s/fulcrum' % (cell,timeP,refControl)
if not os.path.exists(work_dir):
	os.mkdir(work_dir)
# gcto=gct.GCT()
# gcto.read(gctfile)

dp = doseClass.DosePlate()
dp.add_from_gct(gctfile)
dp.examine_doses_tested()
dp.match_template()

#plot whole correlation profile vs. just significant probes
for pert in dp.perts_at_dose:
	templ = 'linear'
	isigCorrs = dp.templateMatchInd[pert][templ]
	corrVec = dp.probe_template_corrs[pert][templ]
	sigCorrs = corrVec[isigCorrs]
	for ind in dp.doseIndDict[pert]:
		doseStr  = dp.gct.get_column_meta('pert_dose')[ind]
		zVec = dp.gct.matrix[:,ind]
		iZsort = np.argsort(zVec)
		zVecSort = zVec[iZsort]
		corrVecSort = corrVec[iZsort]
		plt.plot(zVecSort,corrVecSort,'go')
		plt.xlim(-10,10)
		plt.ylim(-1,1)
		plt.xlabel('z-score')
		plt.ylabel('template corr')
		plt.title(templ)
		outFile = os.path.join(work_dir, pert + '_' + doseStr + '_z_template_corr.png')
		plt.savefig(outFile)
		plt.close()

#plot just signficicant corr markers
for pert in dp.perts_at_dose:
	templ = 'linear'
	isigCorrs = dp.templateMatchInd[pert][templ]
	if len(isigCorrs) > 1:
		corrVec = dp.probe_template_corrs[pert][templ]
		sigCorrs = corrVec[isigCorrs]
		markUp = np.zeros_like(sigCorrs)
		markUp[sigCorrs>0] = sigCorrs[(sigCorrs>0)]
		markDn = np.zeros_like(sigCorrs)
		markDn[sigCorrs<0] = sigCorrs[(sigCorrs<0)]
		markUp = np.multiply(markUp,200)
		markDn = np.multiply(np.abs(markDn),200)
		for count, ind in enumerate(dp.doseIndDict[pert]):
			doseStr  = dp.gct.get_column_meta('pert_dose')[ind]
			zVec = dp.gct.matrix[isigCorrs,ind]
			constant = np.repeat(count +1,len(zVec))
			plt.scatter(zVec,constant,s=markUp,alpha=0.2,c='b')
			plt.scatter(zVec,constant,s=markDn,alpha=0.2,c='r')
		plt.ylim(0,count+2)
		plt.title(pert)
		plt.ylabel('dose um')
		plt.xlabel('z-score')
		doseList = [dp.gct.get_column_meta('pert_dose')[x] for x in dp.doseIndDict[pert]]
		plt.yticks(range(1, len(dp.doseIndDict[pert]) + 1), doseList, rotation = 0)
		outFile = os.path.join(work_dir, pert + '_' + doseStr + '_1D.png')
		plt.savefig(outFile)
		plt.close()

#plot all probes
for pert in dp.perts_at_dose:
	templ = 'linear'
	corrVec = dp.probe_template_corrs[pert][templ]
	# sigCorrs = corrVec[isigCorrs]
	markUp = np.zeros_like(corrVec)
	markUp[corrVec>0] = corrVec[(corrVec>0)]
	markDn = np.zeros_like(corrVec)
	markDn[corrVec<0] = corrVec[(corrVec<0)]
	markUp = np.multiply(markUp,200)
	markDn = np.multiply(np.abs(markDn),200)
	for count, ind in enumerate(dp.doseIndDict[pert]):
		doseStr  = dp.gct.get_column_meta('pert_dose')[ind]
		zVec = dp.gct.matrix[:,ind]
		constant = np.repeat(count +1,len(zVec))
		plt.scatter(zVec,constant,s=markUp,alpha=0.2,c='b')
		plt.scatter(zVec,constant,s=markDn,alpha=0.2,c='r')
	plt.ylim(0,count+2)
	plt.title(pert)
	plt.ylabel('dose um')
	plt.xlabel('z-score')
	doseList = [dp.gct.get_column_meta('pert_dose')[x] for x in dp.doseIndDict[pert]]
	plt.yticks(range(1, len(dp.doseIndDict[pert]) + 1), doseList, rotation = 0)
	outFile = os.path.join(work_dir, pert + '_' + doseStr + '_1D_all_probes.png')
	plt.savefig(outFile)
	plt.close()

###repeat with leave one out
#plot all probes
doses = dp.gct.get_column_meta('pert_dose')
perts = dp.gct.get_column_meta('pert_id')
rids = dp.gct.get_rids()
for pert in dp.perts_at_dose:
	cid_inds = [i for i,x in enumerate(perts) if pert in x]
	pert_doses = [float(doses[x]) for x in cid_inds]
	tmp_tup = zip(pert_doses,cid_inds)
	tmp_tup.sort()
	pert_doses,cid_inds = zip(*tmp_tup)
	#loop through each observation, leave on out and repeat calcuations
	for count,outInd in enumerate(cid_inds):
		cidsLessOne = list(cid_inds)
		cidsLessOne.remove(outInd)
		dosesLessOne = [float(doses[x]) for x in cidsLessOne]
		pert_data = dp.gct.matrix[:,cidsLessOne]
		templ = 'linear'
		if templ == 'linear':
			template_curve = np.array(dosesLessOne)
		# corrVec = dp.probe_template_corrs[pert][templ]
		cc_list = [stats.pearsonr(pert_data[x,:],template_curve) for x in range(len(rids))]
		rho_vec = [cc_list[x][0] for x in range(len(rids))]
		corrVec = np.array(rho_vec)
		# sigCorrs = corrVec[isigCorrs]
		markUp = np.zeros_like(corrVec)
		markUp[corrVec>0] = corrVec[(corrVec>0)]
		markDn = np.zeros_like(corrVec)
		markDn[corrVec<0] = corrVec[(corrVec<0)]
		markUp = np.multiply(markUp,200)
		markDn = np.multiply(np.abs(markDn),200)
		doseStr  = dp.gct.get_column_meta('pert_dose')[outInd]
		zVec = dp.gct.matrix[:,outInd] #plot correlation data with z-scores which were excluded
		constant = np.repeat(count +1,len(zVec))
		plt.scatter(zVec,constant,s=markUp,alpha=0.1,c='b')
		plt.scatter(zVec,constant,s=markDn,alpha=0.1,c='r')
	plt.ylim(0,count+2)
	plt.title(pert)
	# plt.xlim((-20,20))
	plt.ylabel('dose um')
	plt.xlabel('z-score')
	doseList = [dp.gct.get_column_meta('pert_dose')[x] for x in dp.doseIndDict[pert]]
	plt.yticks(range(1, len(dp.doseIndDict[pert]) + 1), doseList, rotation = 0)
	outFile = os.path.join(work_dir, pert + '_' + doseStr + '_leave_one_out_all_probes.png')
	plt.savefig(outFile)
	plt.close()


#test color vector
colorVec = np.zeros_like(corrVec)
colorVec[corrVec>0] = .1
colorVec[corrVec<0] = .9
markVec = np.abs(corrVec)
markVec = np.multiply(markVec,200)
plt.scatter(zVec,constant,s=markVec,alpha=0.2,c=colorVec,cmap=pylab.cm.RdBu)
plt.xlim((-5,5))

plt.scatter(zVec,constant,s=markUp,alpha=0.3,c='b')
plt.scatter(zVec,constant,s=markDn,alpha=0.3,c='r')
plt.xlim((-5,5))
plt.ylim(0,count+2)
plt.title(pert)
# plt.xlim((-20,20))
plt.ylabel('dose um')
plt.xlabel('z-score')
doseList = [dp.gct.get_column_meta('pert_dose')[x] for x in dp.doseIndDict[pert]]
plt.yticks(range(1, len(dp.doseIndDict[pert]) + 1), doseList, rotation = 0)



### identify significant probes --> template corr | z-score
#leave one out cross validation?

#how can we identify probes which are significantly modulated by a perturbation?
#simply saying z>2 is not rigerous
#does looking at z combined with dose outperform either individually?
