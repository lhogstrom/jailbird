#! /usr/bin/env python
'''
analyze the DOSBIO plates - make sc plots etc
'''

import os
import cmap.io.gct as gct

plate = 'DOSBIO001'
cellLine = 'MCF7'
refControl ='pc'
timeP = '24H'
dataDir = '/xchip/cogs/projects/DOS/data'
gctfile = glob.glob(dataDir +'/%s/%s_%s_%s/by_pert_id_pert_dose/%s_%s_%s_COMPZ.MODZ_SCORE_LM_*x978.gctx' % (refControl,plate,cellLine,timeP,plate,cellLine,timeP))
gctfile = gctfile[0]
# gctfile = '/xchip/obelix/pod/brew/vc/PRISM001_A375_24H/by_pert_id_pert_dose/PRISM001_A375_24H_COMPZ.MODZ_SCORE_LM_n58x978.gctx'
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/%s_%s_%s/fulcrum' % (cell,timeP,refControl)
if not os.path.exists(work_dir):
	os.mkdir(work_dir)

gcto=gct.GCT()
gcto.read(gctfile)

