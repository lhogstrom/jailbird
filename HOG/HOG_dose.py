#!/bin/py
'''
load in HOG plates to examine with dose tools
'''
import numpy
import scipy.stats as stats
import cmap.analytics.sc as sc
import cmap.io.gct as gct
import cmap.analytics.dose 
import cmap.analytics.dose as doseClass
import cmap.util.progress as progress
import glob
import os 

#load in HOG plates dir
data_dir = '/xchip/cogs/projects/HOG/data/brew_by_pert_id_pert_dose'

# cellLst = ['PC3', 'MCF7']
# timeLst = ['6H', '24H']
# platePrefix = 'HOG001'
# cell = 'MCF7'
# timeP = '24H'
# refControl = 'pc' #use pc vs vc controled data
# gctfile = glob.glob('%s/%s/%s_%s_%s/by_pert_id_pert_dose/%s_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (data_dir,refControl,platePrefix,cell,timeP,platePrefix,cell,timeP))
# gctfile = gctfile[0]

### use merged gct
gctfile = '/xchip/cogs/projects/HOG/data/brew_by_pert_id_pert_dose/pc/HOG_merged_MCF7_24H/merged_COMPZ.MODZ_SCORE_LM_n862x978.gctx'

work_dir = '/xchip/cogs/projects/HOG/scratch/%s_%s_%s' % (cell,timeP,refControl)
if not os.path.exists(work_dir):
	os.mkdir(work_dir)
db = gct.GCT() #make a gct object
db.read(gctfile)

# dp = doseClass.DosePlate()
# dp.add_from_gct(gctfile)
# dp.examine_doses_tested()
# dp.match_template()
# dp.permutation_template()
# dp.find_modulated(2,4)

pIDs = db.get_column_meta('pert_id')
pDescs = db.get_column_meta('pert_desc')
pIDset = set(pIDs)
uniqID = list(pIDset)
uniqDescs = []
for pID in uniqID:
	iPert = [i for i,x in enumerate(pIDs) if x == pID]
	pDesc = pDescs[iPert[0]]
	uniqDescs.append(pDesc)

for x in uniqDescs:
	print x

