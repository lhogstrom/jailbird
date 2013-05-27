#!/bin/py
'''
load in HOG plates to examine with dose tools
'''
import scipy.stats as stats
import cmap.analytics.sc as sc
import cmap.io.gct as gct
import cmap.analytics.dose 
import cmap.analytics.dose as doseClass
import cmap.util.progress as progress
import matplotlib.pyplot as plt
import glob
import os 
import numpy as np

work_dir = '/xchip/cogs/projects/HOG/scratch/query_by_dose'
rsltF = '/xchip/cogs/projects/HOG/dose_plate_output-by_pert_id_pert_dose/MCF7/24H/pc/apr25/dose_plate_tool.1366908904964/apr26/my_analysis.sig_query_tool.201304261048123/result_WTCS.LM.COMBINED_n286x206076.gctx'
rslt = gct.GCT()
rslt.read_gctx_matrix(rsltF)
rslt.read_gctx_col_meta(rsltF)
# rslt.read_gctx_row_meta(rsltF)
rids = rslt.get_gctx_rid(rsltF)
cids = rslt.get_cids()

col1 = rslt.matrix[:,-1]
plt.hist(col1,40)

### do bioactive compounds have a more significant affogato profile than DMSOs?
cmpds = [x[:22] for x in cids]
plate_name = [x.split('_')[0] for x in rids]
cell_name = [x.split('_')[1] for x in rids]
pert_name = [x.split(':')[1] for x in rids]
# pert_name_short = [x.split(':')[1][:13] for x in rids]
iMCF7 = [i for i,x in enumerate(cell_name) if x == 'MCF7']

# find self connections
SelfConn = {}
for pert in cmpds:
	SelfConn[pert] = {}
	pert_short = pert[:13]
	ind_cmpd = [i for i,x in enumerate(cids) if pert in x]
	iSelf = [i for i,x in enumerate(pert_name) if pert_short in x]
	for i,ind in enumerate(ind_cmpd):
		colName = cids[ind]
		doseName = float(colName.split('_')[1][:-2])
		col1 = rslt.matrix[:,ind]
		iSort = np.argsort(col1)[::-1]
		iRank = np.argsort(iSort)
		SelfConn[pert][doseName] = iRank[iSelf] 
# graph self connections
for pert in cmpds:
	sortKeys = SelfConn[pert].keys()
	sortKeys.sort()
	sKeysStr = [str(x) for x in sortKeys]
	for i,query in enumerate(sortKeys):
		selfRank = SelfConn[pert][query]
		yVals = np.repeat(i+1,len(selfRank))
		plt.scatter(selfRank,yVals)
	plt.xlim((0, rslt.matrix.shape[0]))
	plt.ylim((0,i+2))
	plt.yticks(range(1, i + 2), sKeysStr, rotation = 0)
	plt.xlabel('affogato query rank')
	plt.ylabel('dose')
	plt.title(pert)
	plt.savefig(os.path.join(work_dir,pert + '_self_connections.png'))
	plt.close()

pert = 'BRD-K87909389-001-01-2'
query = 'BRD-K87909389-001-01-2_3.33UM_HOG003_MCF7_24H'
selfRank = SelfConn[pert][query]
#1D plot of self connections


n_top = 2060
for pert in cmpds:
	# pert = 'BRD-K68202742-001-10-8'
	ind_cmpd = [i for i,x in enumerate(cids) if pert in x]
	topMtrx = np.zeros((n_top,len(ind_cmpd)))
	for i,ind in enumerate(ind_cmpd):
		col1 = rslt.matrix[:,ind]
		iSort = np.argsort(col1)[::-1]
		sortCol = col1[iSort]
		topLst = sortCol[:n_top]
		topMtrx[:,i] = topLst
	#box wisker plot of top connection strength
	plt.boxplot(topMtrx)
	plt.savefig(os.path.join(work_dir,pert + '_top_1Perc_connection_strength.png'))
	plt.close()




#histogram of three doses for a perturbation
# f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
# ax1.hist(rslt.matrix[:,ind_cmpd[0]],40)
# ax1.set_title('Sharing both axes')
# ax2.hist(rslt.matrix[:,ind_cmpd[4]],40)
# ax3.hist(rslt.matrix[:,ind_cmpd[-1]],40)


