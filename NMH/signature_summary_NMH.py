#!/bin/py
'''
create a summary for the NMH signatures
'''
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
import re
import scipy.stats
import pylab
import pandas

def _all_indices(value, qlist):
	'''
	input: 1) string and 2) list - return all indices where the input string matches the item in the list
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

#set 1
# cellLst = ['FIBRNPC', 'NEU', 'NPC']
# timeLst = ['6H', '24H']
# set2
cellLst = ['NEU.KCL']
timeLst = ['24H.4H','6H.4H']


### make SC plots and ss and cc summary table
for cell in cellLst:
	for tim in timeLst:
		### make SC plots
		cellLine = cell
		timeP = tim
		refControl = 'pc' #use pc vs vc controled data
		#use file colapsed by pert_id:
		# gctfile = glob.glob('/xchip/obelix/pod/brew/%s/NMH00*_%s_%s/by_pert_id/NMH00*_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		#use file colapsed by rna_well
		gctfile = glob.glob('/xchip/obelix/pod/brew/%s/NMH00*_%s_%s/by_rna_well/NMH00*_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		gctfile.sort() #sort by replicate number
		for set1 in range(len(gctfile)):
			# set1 = 0
			gctfile1 = gctfile[set1] #choose first replicates
			# work_dir = '/xchip/cogs/projects/NMH/signature_summary/%s_%s_%s' % (cell,timeP,refControl)
			work_dir = '/xchip/cogs/projects/NMH/signature_summary'
			if not os.path.exists(work_dir):
				os.mkdir(work_dir)
			db = gct.GCT() #make a gct object
			db.read(gctfile1)
			### make SC plots
			sco = sc.SC()
			sco.add_sc_from_gctx_meta(gctfile1, verbose=False)
			sco.set_thresh_by_specificity(0.8)
			outf =work_dir + '/'+ cell + '_' + tim + '_set' + str(set1+1) + '_SC.png'
			title1 = cell + '_' + tim + '_set' + str(set1+1)
			# sco.plot(out=outf,title=title1)
			fquads = work_dir + '/' + cell + '_' + tim + '_set' + str(set1+1)
			sco.write_all_quad_perts_to_file(fquads)
			# plt.close()
			# ### create table - perturbations ranked by SS/CC
			for metric in ['distil_ss', 'distil_cc_q75']:	
				testVals = numpy.array(db.get_column_meta(metric))
				iSort = numpy.argsort(testVals)[::-1] #sort summary statistic in decending order
				fname = cell + '_' + tim + '_set' + str(set1+1) + '_' + metric + '_rank.txt'
				with open(os.path.join(work_dir,fname),'w') as f:
					headers = ['pert_id','pert_desc','distil_ss',
							   'distil_cc_q75']
					f.write('\t'.join(headers) + '\n')
					for i,pert in enumerate(db.get_column_meta('pert_desc')):
						#prog.update('making SC summary',i,num_perts)
						iRankPosit = iSort[i] #list by ss rank
						data = [db.get_column_meta('pert_id')[iRankPosit],db.get_column_meta('pert_desc')[iRankPosit],db.get_column_meta('distil_ss')[iRankPosit],
							db.get_column_meta('distil_cc_q75')[iRankPosit]]
						f.write('\t'.join(data) + '\n')
			pertIDs = db.get_column_meta('pert_id')
			iDMSO = _all_indices('DMSO',pertIDs)
			#make histogram of ss
			ssStr = db.get_column_meta('distil_ss')
			ss = [float(x) for x in ssStr]
			ssDMSO = [ss[i] for i in iDMSO]
			h1 = plt.hist(ss,30,color='b',range=[0,20],label=['compound'])
			h2 = plt.hist(ssDMSO,30,color='r',range=[0,20],label='DMSO')
			plt.legend()
			plt.xlabel('signature strength')
			plt.ylabel('perterbations per plate')
			plt.title(title1)
			outf =work_dir + '/'+ cell + '_' + tim + '_set' + str(set1+1) + '_ss_hist.png'
			plt.savefig(outf, bbox_inches=0)
			plt.close()
			#make hist of cc
			ccStr = db.get_column_meta('distil_cc_q75')
			cc = [float(x) for x in ccStr]
			ccDMSO = [cc[i] for i in iDMSO]
			plt.hist(cc,30,color='b',range=[-.5,1],label=['compound'])
			plt.hist(ccDMSO,30,color='r',range=[-.5,1],label='DMSO')
			plt.legend()
			plt.xlabel('replicate correlation, q75')
			plt.ylabel('perterbations per plate')
			plt.title(title1)
			outf =work_dir + '/'+ cell + '_' + tim + '_set' + str(set1+1) + '_cc_hist.png'
			plt.savefig(outf, bbox_inches=0)
			plt.close()
			#remove null values
			i666 = [i for i in range(len(cc)) if cc[i] == -666]
			[ss.pop(i) for i in i666]
			[cc.pop(i) for i in i666]
			plt.scatter(ss,cc)
			plt.plot(ss,cc)
			


ssDict = {} #this will contain the signature strength data for each compound on each plate
ccDict = {}
for cell in cellLst:
	for tim in timeLst:
		### make SC plots
		cellLine = cell
		timeP = tim
		refControl = 'pc' #use pc vs vc controled data
		#use file colapsed by pert_id:
		# gctfile = glob.glob('/xchip/obelix/pod/brew/%s/NMH00*_%s_%s/by_pert_id/NMH00*_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		#use file colapsed by rna_well
		gctfile = glob.glob('/xchip/obelix/pod/brew/%s/NMH00*_%s_%s/by_rna_well/NMH00*_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		gctfile.sort() #sort by replicate number
		for set1 in range(len(gctfile)):
			# set1 = 0
			gctfile1 = gctfile[set1] #choose first replicates
			# work_dir = '/xchip/cogs/projects/NMH/signature_summary/%s_%s_%s' % (cell,timeP,refControl)
			work_dir = '/xchip/cogs/projects/NMH/signature_summary' 
			if not os.path.exists(work_dir):
				os.mkdir(work_dir)
			db = gct.GCT() #make a gct object
			db.read(gctfile1)
			#add to ssDict
			plateNameLong = db.get_column_meta('rna_plate')[0]
			plateName = plateNameLong.split('|')[0]
			plateName = plateName.replace('_X1','') #remove replicate number from plate name
			plateName = plateName.replace('_X2','')
			plateName = plateName.replace('_X2','')
			#make nested dictionary
			ssDictSmall = {}
			for i,pert in enumerate(db.get_column_meta('pert_id')):
				ssDictSmall[pert] = float(db.get_column_meta('distil_ss')[i])
			ssDict[plateName] = ssDictSmall
ssFrame = pandas.DataFrame(ssDict)

##access data in ssFrame
for x in ssFrame.index:
	CmpdLst = list(ssFrame.ix[x])


#boxplot - ss by cell line
cellMtrx = numpy.zeros(shape=(len(ssFrame.columns),len(ssFrame.index)))
for i,x in enumerate(ssFrame.columns):
	colLst = list(ssFrame[x])
	cellMtrx[i,:] = colLst
### signature stength of ss in the different plates
data = []
nameLst = []
for j,x in enumerate(ssFrame.columns):
	colLst = [cellMtrx[j,i] for i in range(cellMtrx.shape[1]) if not numpy.isnan(cellMtrx[j,i])] #remove nans
	# plt.boxplot(colLst)
	data.append(colLst)
	nameLst.append(x)
fig = plt.figure()
ax = fig.add_subplot(111)
xtickNames = plt.setp(ax, xticklabels=nameLst)
plt.setp(xtickNames, rotation=45, fontsize=8)
plt.boxplot(data)
plt.show()

### loop through each compound, rank order by mean ss - how consistant is this?
for j,x in enumerate(ssFrame.index):
	rowLst = list(ssFrame.ix[x])
	colLst = [cellMtrx[j,i] for i in range(cellMtrx.shape[1]) if not numpy.isnan(cellMtrx[j,i])] #remove nans


# #load in data by well:
# # fname2 = '/xchip/obelix/pod/brew/pc/NMH001_NEU.KCL_6H.4H/NMH001_NEU.KCL_6H.4H_ZSPCQNORM_n750x978.gctx'
# # fname2 = '/xchip/obelix/pod/brew/pc/NMH001_NEU.KCL_6H.4H/by_pert_id/NMH001_NEU.KCL_6H.4H_COMPZ.MODZ_SCORE_LM_n336x978.gctx'
# fname2 = '/xchip/obelix/pod/brew/pc/NMH001_NEU.KCL_6H.4H/by_rna_well/NMH001_NEU.KCL_6H.4H_COMPZ.MODZ_SCORE_LM_n380x978.gctx'
# db2 = gct.GCT()
# db2.read(fname2)
# pertDescs = db2.get_column_meta('pert_desc')
# iDMSOs = _all_indices('DMSO',pertDescs)
# #can i make an sc plot from the well data?
# sco = sc.SC()
# sco.add_sc_from_gctx_meta(fname2, verbose=False)
# sco.set_thresh_by_specificity(0.8)
# sco.plot()

# pertTypes = db2.get_column_meta('pert_type')
# iCtlV = _all_indices('ctl_vehicle',pertTypes)
# iCtlUT = _all_indices('ctl_untrt',pertTypes)




### junk
# if ssDict.has_key('plate'):
# 	ssDict['plate'].append(plateName)
# 	ccDict['plate'].append(plateName)
# else:
# 	ssDict['plate'] = [plateName]
# 	ccDict['plate'] = [plateName]
# for metric in ['distil_ss', 'distil_cc_q75']:
# 	testVals = numpy.array(db.get_column_meta(metric))
# 	# iSort = numpy.argsort(testVals)[::-1] #sort summary statistic in decending order
# 	for i,pert in enumerate(db.get_column_meta('pert_id')):
# 		# iRankPosit = iSort[i] #list by ss rank
# 		if metric == 'distil_ss':
# 			if ssDict.has_key(pert):
# 				ssDict[pert].append(db.get_column_meta('distil_ss')[i])
# 			else:
# 				ssDict[pert] = [db.get_column_meta('distil_ss')[i]]
# 		if metric == 'distil_cc_q75':
# 			if ccDict.has_key(pert):
# 				ccDict[pert].append(db.get_column_meta('distil_cc_q75')[i])
# 			else:
# 				ccDict[pert] = [db.get_column_meta('distil_cc_q75')[i]]
# ssDict1 = {}
# ssDict1[plateName] = [db.get_column_meta('distil_ss')]
# for i,pert in enumerate(db.get_column_meta('pert_id')):
# 	ssDict1[pert] = float(db.get_column_meta('distil_ss')[i])
# #make data frame of ss and cc data
# n_perts = len(db.get_column_meta('pert_id'))

# ssFrame[plateName] = ssSeries
# ssFrame = pandas.DataFrame(ssDict1,index=ssDict1.keys(),columns=[plateName])
# ssFrame = pandas.DataFrame(columns=[plateName],)

# ssDict1 = {}
# ssDict1[plateName] = db.get_column_meta('distil_ss')
# ssFrame = pandas.DataFrame(ssDict1,index=db.get_column_meta('pert_id'))
# #add another set of data
# ssSeries = pandas.Series(ssDict1)
# ssFrame['plate2'] = ssSeries


#heatmaps
#qq plots

#test quads out
db.read(gctfile1)
### make SC plots
sco = sc.SC()
sco.add_sc_from_gctx_meta(gctfile1, verbose=False)
sco.set_thresh_by_specificity(0.8)
outf =work_dir + '/'+ cell + '_' + tim + '_set' + str(set1+1) + '_SC.png'
title1 = cell + '_' + tim + '_set' + str(set1+1)
# sco.plot(out=outf,title=title1)
fquads = work_dir + '/' + cell + '_' + tim + '_set' + str(set1+1)
sco.write_all_quad_perts_to_file(fquads)

work_dir = '/xchip/cogs/projects/NMH/signature_summary'
gctfile1 = '/xchip/obelix/pod/brew/pc/NMH002_NEU.KCL_6H.4H/by_rna_well/NMH002_NEU.KCL_6H.4H_COMPZ.MODZ_SCORE_LM_n380x978.gctx'
sco = sc.SC()
sco.add_sc_from_gctx_meta(gctfile1, verbose=False)
sco.set_thresh_by_specificity(0.8)
fquads = work_dir + '/' + cell + '_' + tim + '_set' + str(set1+1)
sco.write_all_quad_perts_to_file(fquads)


