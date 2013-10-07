#!/bin/py
'''
pull compounds in the top right of the SC plot and cluster them
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


### make top right lists
for cell in cellLst:
	for tim in timeLst:
		### make SC plots
		cellLine = cell
		timeP = tim
		refControl = 'vc' #use pc vs vc controled data
		#use file colapsed by pert_id:
		# gctfile = glob.glob('/xchip/obelix/pod/brew/%s/NMH00*_%s_%s/by_pert_id/NMH00*_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		#use file colapsed by rna_well
		gctfile = glob.glob('/xchip/obelix/pod/brew/%s/NMH00*_%s_%s/by_rna_well/NMH00*_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
		gctfile.sort() #sort by replicate number
		for set1 in range(len(gctfile)):
			# set1 = 0
			gctfile1 = gctfile[set1] #choose first replicates
			# work_dir = '/xchip/cogs/projects/NMH/signature_summary/%s_%s_%s' % (cell,timeP,refControl)
			work_dir = '/xchip/cogs/projects/NMH/top_right_cluster'
			if not os.path.exists(work_dir):
				os.mkdir(work_dir)
			db = gct.GCT() #make a gct object
			db.read(gctfile1)
			### make SC plots
			sco = sc.SC()
			sco.add_sc_from_gctx_meta(gctfile1, verbose=False)
			sco.set_thresh_by_specificity(0.8)
			# outf =work_dir + '/'+ cell + '_' + tim + '_set' + str(set1+1) + '_SC.png'
			# title1 = cell + '_' + tim + '_set' + str(set1+1)
			# sco.plot(out=outf,title=title1)
			fquads = work_dir + '/' + 'NMH00' + str(set1+1) + '_' + cell + '_' + tim
			sco.write_all_quad_ids_to_file(fquads)



