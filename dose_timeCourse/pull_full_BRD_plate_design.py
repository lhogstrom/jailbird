#!/bin/py
import numpy
import os
import cmap.util.mongo_utils as mutil
import cmap.io.gct as gct
import cmap.analytics.sc as sc
import matplotlib.pyplot as plt

#set output dir
work_dir = '/xchip/cogs/projects/ASG_dose_time/plate_design/short_brd'
fname = os.path.join(work_dir,'parp_list.txt')
tstPanLst = open(fname).read().splitlines()
pIDsFull = []
for pID in tstPanLst:
	cpcLst = []
	CM = mutil.CMapMongo()
	brewPrfx = 'CPC006'
	#query only CPC plates
	# cpcLst = CM.find({'brew_prefix':{'$regex':brewPrfx},'pert_id':{'$regex':pID}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
	# cpcLst = CM.find({'brew_prefix':{'$regex':brewPrfx},'pert_id':{'$regex':pID}},{'pert_id':True},limit=1)
	# query all of affogato
	cpcLst = CM.find({'pert_id':{'$regex':pID}},{'pert_id':True},limit=1)
	if cpcLst:
		pIDsFull.append(cpcLst[0])
	else:
		print 'no match for' + pID
	# for cpc in cpcLst:
	# 	if len(cpc) > 21:
	# 		pIDsFull.append(cpc)
	# 		if 
	# 		continue

### get brds for pert_descs
fname = os.path.join(work_dir,'parp_common_names.txt')
tstPanLst = open(fname).read().splitlines()
pIDsFull = []
for pID in tstPanLst:
	cpcLst = []
	CM = mutil.CMapMongo()
	cpcLst = CM.find({'pert_desc':{'$regex':pID}},{'pert_id':True},limit=1)
	if cpcLst:
		pIDsFull.append(cpcLst[0])
	else:
		print 'no match for' + pID

for x in pIDsFull:
	print x




# query for individual brds
pID = 'BRD-A75409952'
CM = mutil.CMapMongo()
brewPrfx = 'CPC006'
#query only CPC plates
# cpcLst = CM.find({'brew_prefix':{'$regex':brewPrfx},'pert_id':{'$regex':pID}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
# cpcLst = CM.find({'brew_prefix':{'$regex':brewPrfx},'pert_id':{'$regex':pID}},{'pert_id':True},limit=1)
# query all of affogato
cpcLst = CM.find({'pert_id':{'$regex':pID}},{'pert_id':True})


# check overlap between groups
fname = os.path.join(work_dir,'pan_list_plate1.txt')
plate1 = open(fname).read().splitlines()
fname2 = os.path.join(work_dir,'pan_list_plate2.txt')
plate2 = open(fname2).read().splitlines()
set1 = set(plate1)
set2 = set(plate2)

for x in set2:
	print x

### get pert_descs that match with each Full BRD
fname = os.path.join(work_dir,'full_brds.txt')
tstPanLst = open(fname).read().splitlines()
pdescDict = {}
for pID in tstPanLst:
	cpcLst = []
	# CM = mutil.CMapMongo(mongo_location='a2y13q1')
	CM = mutil.CMapMongo()
	cpcLst = CM.find({'pert_id':pID},{'pert_desc':True})
	if cpcLst:
		pdescDict[pID] = cpcLst
		cpcLst = []
	else:
		pdescDict[pID] = []
		print 'no match for' + pID
uniquePdesc = {}
for pID in pdescDict:
	uniquePdesc[pID] = set(pdescDict[pID])
fout = os.path.join(work_dir,'id_descs.txt')
with open(fout,'w') as f:
	for pID in tstPanLst:
		f.write(pID + '\t')
		for desc in uniquePdesc[pID]:
			f.write(desc + '\t')
		f.write('\n')

### get pert_descs that match with each SHORT BRD
### FOR SHORT BRD
fname = os.path.join(work_dir,'full_brds.txt')
tstPanLst = open(fname).read().splitlines()
pdescDict = {}
for pID in tstPanLst:
	pID = pID[:13] 
	cpcLst = []
	# CM = mutil.CMapMongo(mongo_location='a2y13q1')
	CM = mutil.CMapMongo()
	cpcLst = CM.find({'pert_id':{'$regex':pID}},{'pert_desc':True})
	if cpcLst:
		pdescDict[pID] = cpcLst
		cpcLst = []
	else:
		pdescDict[pID] = []
		print 'no match for' + pID
uniquePdesc = {}
for pID in pdescDict:
	pID = pID[:13]
	uniquePdesc[pID] = set(pdescDict[pID])
fout = os.path.join(work_dir,'id_descs.txt')
with open(fout,'w') as f:
	for pID in tstPanLst:
		pID = pID[:13]
		f.write(pID + '\t')
		for desc in uniquePdesc[pID]:
			if desc == '-666.0':
				continue
			else:
				f.write(desc + '\t')
		f.write('\n')

