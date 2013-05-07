#! /usr/bin/env python
'''
work out what individual dose pages will look like
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
import cmap.analytics.fitting as fitting
import glob
import subprocess
import re

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
# work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/ec50_test_inf/%s_%s_%s' % (cell,timeP,refControl)
# if not os.path.exists(work_dir):
# 	os.mkdir(work_dir)
gcto = gct.GCT() #make a gct object
gcto.read(gctfile)

work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/A375/6H/dose_plate_tool/apr03/dose_plate_tool.1365019286965'
# gcto = gct.GCT(args.res)
# gcto.read()
do = DosePlate()
# do.add_from_gct(args.res)
do.add_from_gct(gctfile)
do.examine_doses_tested()
cids = gcto.get_cids()
pert_descs = gcto.get_column_meta('pert_desc')
perts = [x.split(':')[1] for x in cids]
pert_desc_dict = dict(zip(perts,pert_descs))
unique_perts = list(set(perts))
unique_perts.sort()
unique_pDescs = list(pert_desc_dict.values())
unique_pDescs.sort()

#Make a nested dictionary with the following structure:
# key = pertID
# values (nested lists):
	# 1) dose concentrations tested for each compound
	# 2) self connection percenatge
	# 3) signature strength 
	# 4) signature strength as percentage of maximum
	# 5) replicate correlation 
	# 6) compound subpage html
	# 7) number of doses tested
topPercDic = {} #dict - key = pert_desc, value = [top query percent for each dose], [ss], [ss as percent of max], [cc], [compound page link]
qSumF = os.path.join(work_dir, 'plate_self_connection.txt')
counter = 0
with open(qSumF,'rt') as f:
	for line in f:
		splt = line.split()
		if counter == 0: #get dose data from first line and skip
			for i,y in enumerate(splt):
				if i >= 1: #don't shorten the name
					splt[i] = y[:5] #shorten to 3 sig-figs
				doseUnits = splt
			counter = 1
			continue
		for i,y in enumerate(splt):
			if i >= 1: #don't shorten the name
				splt[i] = y[:5] #shorten to 3 sig-figs
		topPercDic[splt[0]] = [doseUnits[1:]] #make first entry the dose units tested
		topPercDic[splt[0]].extend([splt[1:]]) #make second entry the self connection
# load in signature strength 
ssF = os.path.join(work_dir, 'sig_strength_summary.txt')
counter = 0
with open(ssF,'rt') as f:
	for line in f:
		if counter == 0: #skip first line
			counter = 1
			continue
		splt = line.split()
		for i,y in enumerate(splt):
			if i >= 1: #don't shorten the name
				splt[i] = y[:5] #shorten to 3 sig-figs
		topPercDic[splt[0]].extend([splt[1:]])
#load in signature strength as percentage of maximum
ssMaxF = os.path.join(work_dir, 'sig_strength_percentOfMax.txt')
counter = 0
with open(ssMaxF,'rt') as f:
	for line in f:
		if counter == 0: #skip first line
			counter = 1
			continue
		splt = line.split()
		for i,y in enumerate(splt):
			if i >= 1: #don't shorten the name
				splt[i] = y[:5] #shorten to 3 sig-figs
		topPercDic[splt[0]].extend([splt[1:]])
#put CC data into value for each compound
ccF = os.path.join(work_dir, 'cc_summary.txt')
counter = 0
with open(ccF,'rt') as f:
	for line in f:
		if counter == 0: #skip first line
			counter = 1
			continue
		splt = line.split()
		for i,y in enumerate(splt):
			if i >= 1: #don't shorten the name
				splt[i] = y[:5] #shorten to 3 sig-figs
		topPercDic[splt[0]].extend([splt[1:]])
# add link to dose dictionary
# for x in unique_pDescs:
for x in topPercDic:
	if x == 'DMSO':
		continue
	else:
		lnk = x + '_detail.html'
		topPercDic[x].append(lnk)
#add number of doses tested to dose dictionary
# for x in unique_pDescs:
for x in topPercDic:
	if x == 'DMSO':
		continue
	else:
		nDoses = len(topPercDic.values()[0][0])
		topPercDic[x].append(nDoses)


# buld an environment for jinja2
cmap_base_dir = '/'.join(os.path.dirname(cmap.__file__).split('/')[0:-1])
env = jinja2.Environment(loader=jinja2.FileSystemLoader(cmap_base_dir + '/templates'))
index_page_template = env.get_template('dose_plate_tool_Template.html')
index_links = [x + '_detail.html' for x in unique_pDescs]
with open(os.path.join(work_dir,'index.html'),'w') as f:
	f.write(index_page_template.render(title='Dose Analysis Results',
										data=topPercDic))	

#make sub-pages
fnames = os.listdir((os.path.join(work_dir))) # list of all files in working directory
# for each unique_pert, make a detail page
dose_response_compound_summary_template = env.get_template('Dose_Response_Compound_Summary_Template.html')
for unique_pert in do.perts_at_dose:
	scImages = []
	for l in fnames: #make list of query rank images 
		matchObj = re.search(unique_pert + '(.*)_SC.png(\.*)', l, re.M|re.I)
		if matchObj:
			scImages.append(matchObj.group())
	qq_images = []
	for l in fnames: #make list of query rank images 
		matchObj = re.search(unique_pert + '(.*)_qq.png(\.*)', l, re.M|re.I)
		if matchObj:
			qq_images.append(matchObj.group())
	heatmap_images = []
	for l in fnames: #make list of query rank images 
		matchObj = re.search(unique_pert + '(.*)_heatmap.png(\.*)', l, re.M|re.I)
		if matchObj:
			heatmap_images.append(matchObj.group())
	with open(os.path.join(work_dir,unique_pert + '_detail.html'),'w') as f:
		f.write(dose_response_compound_summary_template.render(
				title=unique_pert,
				scImages=scImages,
				heatmap_images=heatmap_images,
				qq_images=qq_images))

