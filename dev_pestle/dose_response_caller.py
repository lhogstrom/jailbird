import cmap.analytics.sc as sc
import cmap.io.gct as gct
import cmap.analytics.dose as doseClass
import cmap.util.progress as progress
import scipy.stats as stats
import numpy as np
import copy
import os
import matplotlib.pyplot as plt
import tables

def dose_response_caller(gct_file,verbose=True):
	'''
	Counts the number of dose-responsive genes based on template matching
	also examines the S-C correlation to determine if a perturbation is dose 
	responsive
	'''
	dp = doseClass.DosePlate()
	dp.add_from_gct(gct_file)
	dp.examine_doses_tested(verbose=verbose)
	dp.match_template(verbose=verbose)
	# dp.permutation_template()
	ssAll = dp.gct.get_column_meta('distil_ss')
	ccAll = dp.gct.get_column_meta('distil_cc_q75')
	#matric what is the dose corr
	nMatchDict = {}
	for pert in dp.perts_at_dose:
		IndsPert = dp.doseIndDict[pert]
		ssVals = [float(ssAll[i]) for i in IndsPert]
		ccVals = [float(ccAll[i]) for i in IndsPert]
		corr, p = stats.pearsonr(ccVals,ssVals)
		if p < .2:
			print pert + ' S-C corr ' + str(corr)
		#find number of probes that match a template pass FDR
		templates = dp.templateMatchInd[pert].keys()
		probeTempMatches = []
		for temp in templates:
			probeTempMatches.extend(dp.templateMatchInd[pert][temp])
		probeMatch = set(probeTempMatches)
		print pert + ' number of template matches ' + str(len(probeMatch))
		### DR caller
		if len(probeMatch) > 10 or p < .1:
			print pert + ' is dose responsive'
		nMatchDict[pert] = len(probeMatch)
	return nMatchDict

# execfile('/xchip/cogs/hogstrom/scripts/dev_pestle/dose_class_obj.py')
#brew by rna_well
# gctfile = '/xchip/obelix/pod/brew_tmp/pc/PRISM001_PC3_6H/by_rna_well/PRISM001_PC3_6H_COMPZ.MODZ_SCORE_LM_n260x978.gctx'
#brew by pert_id_pert_dose
gctfile = '/xchip/cogs/data/brew/a2y13q1/PRISM001_PC3_24H/by_pert_id_pert_dose/PRISM001_PC3_24H_COMPZ.MODZ_SCORE_LM_n60x978.gctx'
dose_response_caller(gctfile)


### permute columns of plate gct files to check caller
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/dose_response_call'
db = gct.GCT()
db.read(gctfile)
randGCT = copy.copy(db)
n_random = 100
nMatchList = {}
for rand in range(n_random):
	### shuffle columns and rows sepratly
	# indShuffle1 = range(db.matrix.shape[1])
	# indShuffle2 = range(db.matrix.shape[0])
	# np.random.shuffle(indShuffle1)
	# np.random.shuffle(indShuffle2)
	# randGCT.matrix = db.matrix[:,indShuffle1] #perumute columns 
	# randGCT.matrix = randGCT.matrix[indShuffle2,:] #perumute rows columns 
	### shuffle all data toghether
	randMtrx = db.matrix.flatten()
	indShuffle = range(len(randMtrx))
	np.random.shuffle(indShuffle)
	randMtrx = randMtrx[indShuffle]
	randGCT.matrix = randMtrx.reshape(db.matrix.shape)
	## jumble ss and cc column annotations ## 
	hdf5Table = tables.openFile(randGCT.src,'a') #open file with append
	anntTable = hdf5Table.listNodes("/0/META/COL")
	# ccTable = [x for x in anntTable[12].iterrows()]
	i_trt_cp = [i for i,x in enumerate(db.get_column_meta('pert_type')) if x == 'trt_cp'] #index of all trt_cp
	iRandCol = i_trt_cp[:]
	np.random.shuffle(iRandCol)
	#scrable values for cc and ss for trt_cps
	for i2,i1 in enumerate(i_trt_cp):
		anntTable[16][i1] = anntTable[16][iRandCol[i2]]
		anntTable[12][i1] = anntTable[12][iRandCol[i2]]
	hdf5Table.close()
	#write matrix
	randF = os.path.join(work_dir,'permuted_data.gctx')
	# randGCT._close_gctx()
	randGCT.src = randF
	randGCT.write(randF)
	nMatch = dose_response_caller(randF,verbose=False)
	# nMatchList.append(nMatch)
	for x in nMatch:
		if nMatch[x] > 100:
			print x + ', ' + nMatch[x] + ' this is huge'
		if nMatchList.has_key(x):
			nMatchList[x].append(nMatch[x])
		else:
			nMatchList[x] = []
			nMatchList[x].append(nMatch[x])

combineMatches = []
for pert in nMatchList.keys():
	combineMatches.extend(nMatchList[pert])


### create empty pytable file
# hdf = tables.openFile('test.h5','w',title='hi pytables')
# hdf.createGroup('/','box')
# data = np.arange(24).reshape((8,3))
# grid = hdf.createArray('/box', 'grid', data)
# #modify grid entries
# grid[0] = [10, 11, 12]	
# hdf.renameNode('/box','folder')
# hdf.copyNode(grid,'/')
# hdf.root._v_attrs.author = "Ivan"

# data2 = hdf.listNodes('/box')[0]
# data2[2] = [33, 33, 33]

# objectTree 
# --> root group
# group = node
# leaf = datasets
# nodes created through methods

### plot heatmap
# pert1 = 'BRD-K01737880'
# ind880 = np.array(dp.doseIndDict[pert1])
# flagged = dp.templateMatchInd[pert1]['linear']
# plt.imshow(randGCT.matrix[flagged,ind880],interpolation='nearest',cmap='RdBu')
# plt.imshow(g[:10],interpolation='nearest')
