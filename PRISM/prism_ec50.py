#!/bin/py
'''
perform analysis ec50 with prism
'''
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import collections
import pandas as pd
import cmap.io.gct as gct
import cmap.util.tool_ops as tool_ops
import cmap.util.progress as progress
import cmap.analytics.fitting as fitting
import cmap.analytics.cluster as cluster
import cmap.util.queue as queue
import cmap.io.plategrp as grp
import cmap.plot.colors as colors
import cmap.io.rnk as rnk
import cmap.analytics.es as es
import glob

#return all indices where the input string matches the item in the list
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

#load in prism plates
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism'

#run with all cell lines/ time points
# cellLst = ['PC3', 'A375', 'MCF7']
# timeLst = ['6H', '24H']
cell = 'PC3'
tim = '6H'
cellLine = cell
timeP = tim
refControl = 'pc' #use pc vs vc controled data
gctfile = glob.glob('/xchip/obelix/pod/brew/%s/PRISM001_%s_%s/by_pert_id_pert_dose/PRISM001_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (refControl,cellLine,timeP,cellLine,timeP))
gctfile = gctfile[0]
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/ec50_test/%s_%s_%s' % (cell,timeP,refControl)
if not os.path.exists(work_dir):
	os.mkdir(work_dir)
db = gct.GCT() #make a gct object
db.read(gctfile)

### def ec50(args, work_dir):
# gcto = gct.GCT(args.res)
# gcto.read()
gcto = db

# grab the cids from the file and mine dose information from them.  Find all of
# the unique perts
cids = gcto.get_gctx_cid(gctfile)
pert_descs = gcto.get_column_meta('pert_desc')
doses = [float(x.split(':')[2]) for x in cids]
perts = [x.split(':')[1] for x in cids]
unique_perts = list(set(perts))

# grab the cid and rid for use below
rids = gcto.get_gctx_rid(gctfile)

input_list = []
for i, unique_pert in enumerate(unique_perts):
    cid_inds = [i for i, x in enumerate(cids) if unique_pert in x]
    pert_doses = [doses[x] for x in cid_inds]
    tmp_tup = zip(pert_doses, cid_inds)
    tmp_tup.sort()
    pert_doses, cid_inds = zip(*tmp_tup)
    tup = (cids,unique_pert,pert_descs,doses,gcto,rids,work_dir)
    input_list.append(tup)

### def _par_plot(input,output):
# for input_tup in iter(input.get, 'STOP'):
for input_tup in input_list:
    # unpack the list of passed variables from the input
    cids = input_tup[0]
    unique_pert = input_tup[1]
    pert_descs = input_tup[2]
    doses = input_tup[3]
    gcto = input_tup[4]
    rids = input_tup[5]
    work_dir = input_tup[6]

    # grab the z-scores and doses for the current pert and sort the pairs
    # by dose. put the cid_inds in the same sorted order
    cid_inds = [i for i,x in enumerate(cids) if unique_pert in x]
    pert_desc = pert_descs[cid_inds[0]] #set pert desc to the first dose
    pert_desc = pert_desc.replace(' ','_')
    pert_doses = [doses[x] for x in cid_inds]
    tmp_tup = zip(pert_doses,cid_inds)
    tmp_tup.sort()
    pert_doses,cid_inds = zip(*tmp_tup)

    print('started {0}'.format(pert_desc))
    if len(pert_doses) > 1:
        # compute the ec50 of all probes for the current pert
        cols = [cids[x] for x in cid_inds]
        pert_df = gcto.frame[cols]
        pos_mod_ind = set(pert_df[pert_df.max(axis=1) > 2].index) #list probes that were above 2 z - at any dose?
        neg_mod_ind = set(pert_df[pert_df.min(axis=1) < -2].index)
        mod_ind = pos_mod_ind.union(neg_mod_ind)
        pert_df = pert_df.reindex(mod_ind)
        pert_data = pert_df.values
        EC50 = fitting.EC50(list(pert_doses),pert_data) #fit data to 
        EC50.compute_ec50(Logarithmic=False)
        #EC50.ec50.shape - number of active probes
        #EC50.ec50_matrix.shape - number of probes x 4?
        #EC50.hill_coeficient - what is the hill coeff
        #EC50.responses - #doses x # of active probes - a list of lists or 2D array of expression responses
        #core 	line: v0 = self._calc_init_params(self.concentrations,response)
        #		def _calc_init_params(self,x,y):
        #				'''
		#		 This generates the min, max, x value at the mid-y value, and Hill
		# 		 coefficient. These values are starting points for the sigmoid fitting.
		# 		 x & y are the points to be fit
		# 		 returns minimum, maximum, ec50 and hill coefficient starting points
		# 		 '''

        low_pass = EC50.ec50 >= pert_doses[0]
        high_pass = EC50.ec50 <= pert_doses[-1]
        passing_inds = low_pass * high_pass #must have EC50 be within range of doses tested
        non_passing_inds = np.invert(passing_inds) 
        rid_inds = np.arange(len(pert_df.index))
        irrational_rid_inds = rid_inds[non_passing_inds]
        rational_rid_inds = rid_inds[passing_inds]
        rational_ec50 = EC50.ec50[passing_inds]
        rational_pert_df = pert_df.reindex([rids[x] for x in rid_inds[passing_inds]])
        rational_pert_data = rational_pert_df.values
        finite_inds = np.where(np.isfinite(rational_pert_data[:,0]))[0]
        rational_ec50 = rational_ec50[finite_inds]

