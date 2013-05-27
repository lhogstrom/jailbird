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

# # cellLst = ['PC3', 'MCF7']
# # timeLst = ['6H', '24H']
# platePrefix = 'HOG001'
# cell = 'MCF7'
# timeP = '24H'
# refControl = 'pc' #use pc vs vc controled data
# gctfile = glob.glob('%s/%s/%s_%s_%s/by_pert_id_pert_dose/%s_%s_%s_COMPZ.MODZ_SCORE_LM_*.gctx' % (data_dir,refControl,platePrefix,cell,timeP,platePrefix,cell,timeP))
# gctfile = gctfile[0]

# work_dir = '/xchip/cogs/projects/HOG/scratch/%s_%s_%s' % (cell,timeP,refControl)
# if not os.path.exists(work_dir):
# 	os.mkdir(work_dir)
# db = gct.GCT() #make a gct object
# db.read(gctfile)

###merged gct 
gctfile = '/xchip/cogs/projects/HOG/data/brew_by_pert_id_pert_dose/pc/HOG_merged_MCF7_24H/merged_COMPZ.MODZ_SCORE_LM_n862x978.gctx'
db = gct.GCT() #make a gct object
db.read(gctfile)

do = doseClass.DosePlate()
do.add_from_gct()

pertIDs = db.get_column_meta('pert_id')
pertDescs = db.get_column_meta('pert_desc')

### parp group
# XAV939	BRD-K12762134-001-02-1
# veliparib	BRD-K87142802-001-02-7
# olaparib	BRD-K02113016
# 3-amino-benzamide	BRD-K08703257-001-06-3
# IWR-1-endo	BRD-K61314889-001-01-0
# NU-1025	BRD-K48692744-001-01-2
# DR2313	BRD-K94920105-001-01-3
# phenanthridinone	BRD-K11163873-001-02-8
# PJ34	BRD-K11853856-003-01-3
#just BRDs
# BRD-K12762134-001-02-1 x
# BRD-K87142802-001-02-7 x
# BRD-K08703257-001-06-3 x 
# BRD-K61314889-001-01-0 x
# BRD-K48692744-001-01-2 x
# BRD-K94920105-001-01-3 x
# BRD-K11163873-001-02-8 x
# BRD-K11853856-003-01-3 x

parpLst = ['BRD-K12762134-001-02-1', 'BRD-K87142802-001-02-7', 'BRD-K08703257-001-06-3' , 'BRD-K61314889-001-01-0' , 'BRD-K48692744-001-01-2' , 'BRD-K94920105-001-01-3' , 'BRD-K11163873-001-02-8' , 'BRD-K11853856-003-01-3']
kinaseLst = ['BRD-K01737880', 'BRD-K59369769-001-05-4', 'BRD-K36740062-001-02-5', 'BRD-K83963101-001-01-0', 'BRD-K29830875', 'BRD-K63923597']

groupName = 'PARP'
groupName = 'KINASE'
group1 = kinaseLst
iGroup = []
#identify pert IDs that match group
for item in group1:
	item = item[:13]
	iMatch = []
	iMatch = [i for i,ID in enumerate(pertIDs) if item in ID]
	if not iMatch:
		print item + ' not found'
	iGroup.extend(iMatch)
#collect group data and annotations
pert_data = db.matrix[:,iGroup]
group_cids = [db.get_cids()[i] for i in iGroup]
group_pIDs = [db.get_column_meta('pert_id')[i] for i in iGroup]
colDict = {}
for chd in db.get_chd():
	meta = db.get_column_meta(chd)
	colDict[chd] = [meta[i] for i in iGroup]
groupStruc = gct.GCT()
groupStruc.build(pert_data,db.get_rids(),group_cids,{},colDict)
groupStructF = os.path.join(work_dir,groupName + '_compound_group_pc.gctx')
groupStruc.write(groupStructF)


### kinase group
# BRD-K01737880
# BRD-K59369769-001-05-4
# BRD-K36740062-001-02-5
# BRD-K83963101-001-01-0
# BRD-K29830875
# CMAP-AZD-1152HQPA

# write pert_ids to a file
# IdFile = work_dir + 'pert_id_list.txt'
# with open(IdFile,'w') as f:
# 	for ID in pertIDs:
# 		f.write(ID + '\n')