#!/bin/py
import numpy
import os
import cmap.util.mongo_utils as mutil
import cmap.io.gct as gct
import cmap.analytics.sc as sc
import matplotlib.pyplot as plt

#set output dir
work_dir = '/xchip/cogs/projects/ASG_dose_time/plate_design/tmp'
#run mongo query
# qStr = 'DMSO'
cellFirstCheck = ['MCF7','PC3']
candidatePertID = []
candidatePertDesc = []
for cell1 in cellFirstCheck:
	ssMin = 10
	CM = mutil.CMapMongo()
	cmpdLst = CM.find({'cell_id':cell1,'pert_type':'trt_cp','distil_ss':{'$gt':ssMin}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True})
	#is there a way to query mongo for instances with a SS > n? --> 'distil_ss':{'$gt':ssMin}

	#order output by cc
	cmpdCC = [x['distil_cc_q75'] for x in cmpdLst]
	cmpdPertDesc = [x['pert_desc'] for x in cmpdLst]
	cmpdPertID = [x['pert_id'] for x in cmpdLst]
	#sort by CC
	iSortCC = numpy.argsort(cmpdCC)[::-1] #indices that sort CC
	pertDescSort = [cmpdPertDesc[i] for i in iSortCC]
	pertIDSort = [cmpdPertID[i] for i in iSortCC]
	CCSort = [cmpdCC[i] for i in iSortCC]

	#for top hits (with high SS and CC), find if they have high ss in other cell lines
	for i,pID in enumerate(pertIDSort):
		pID = pID[:13]
		cmpd = pertDescSort[i]
		# find all instances of that compound in affogato
		mtchLst = CM.find({'pert_id':{'$regex':pID}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
		ssLst = [x['distil_ss'] for x in mtchLst]
		sigLst = [x['sig_id'] for x in mtchLst]
		cellLst = [x['cell_id'] for x in mtchLst]
		unqCell = set(cellLst)
		nCells = len(unqCell)
		n_sigs = len(sigLst)
		if n_sigs <= 20: #skip if there is less than 10 signatures in affogato 
			continue
		if numpy.mean(ssLst) <= 8: #only examine compounds whos average ss passes a threshold
			continue
		#make S-C plot
		sco = sc.SC()
		sco.add_sc_from_mongo_sig_ids(sigLst)
		# sco.set_thresh_by_specificity(0.8)
		sco.plot(title=(str(n_sigs) + ' signatures of ' + str(cmpd) + 'in affogato, ' + str(nCells) + 'cell lines'),pos_con=['None'],out=os.path.join(work_dir,'_'.join([pID,'SC.png'])))
		candidatePertID.append(pID)
		candidatePertDesc.append(cmpd)
	
#what are the poscons?
CM = mutil.CMapMongo()
cmpdLst = CM.find({'pert_id':{'$regex':qStr},'cell_id':cell1},{'sig_id':True})

### analyze CPC006 plates
work_dir = '/xchip/cogs/projects/ASG_dose_time/plate_design/CPC006/bioactive_tool_list'
# work_dir = '/xchip/cogs/projects/ASG_dose_time/plate_design/CPC006/all_affogato'
CM = mutil.CMapMongo()
brewPrfx = 'CPC006'
cpcLst = CM.find({'brew_prefix':{'$regex':brewPrfx}},{'sig_id':True,'pert_id':True})
sigLst = [x['sig_id'] for x in cpcLst]
pertLst = [x['pert_id'] for x in cpcLst]
uniqueCPCcmpds = set(pertLst)
sco = sc.SC()
sco.add_sc_from_mongo_sig_ids(cpcLst)
#read in output from bioactivity_tool
# fname = '/xchip/cogs/projects/bioactivity/cfwork/CPC006_bioactivity/bioactive_table.txt'
# fname2 = '/xchip/cogs/projects/bioactivity/cfwork/CPC006_bioactivity/test_pan_list.txt'
fname2 = '/xchip/cogs/projects/ASG_dose_time/plate_design/short_brd/pan_list.txt'
tstPanLst = open(fname2).read().splitlines()
tstPanLst = [x.split(':')[0] for x in tstPanLst]
work_dir = '/xchip/cogs/projects/ASG_dose_time/plate_design/pan_SC_plots'

pDescs = []
for pID in tstPanLst:
	CM = mutil.CMapMongo()
	brewPrfx = 'CPC006'
	#query only CPC plates
	cpcLst = CM.find({'brew_prefix':{'$regex':brewPrfx},'pert_id':{'$regex':pID}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
	#query all of db
	# cpcLst = CM.find({'pert_id':{'$regex':pID}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
	ssLst = [x['distil_ss'] for x in cpcLst]
	sigLst = [x['sig_id'] for x in cpcLst]
	cellLst = [x['cell_id'] for x in cpcLst]
	pDescLst = [x['pert_desc'] for x in cpcLst]
	pDescs.append(pDescLst[0])
	pIDLst = [x['pert_id'] for x in cpcLst]
	unqCell = set(cellLst)
	nCells = len(unqCell)
	n_sigs = len(sigLst)
	sco = sc.SC()
	sco.add_sc_from_mongo_sig_ids(sigLst)
	# sco.set_thresh_by_specificity(0.8)
	sco.plot(title=(str(n_sigs) + ' signatures of ' + str(pID) + ' in affogato, ' + str(nCells) + 'cell lines'),pos_con=['None'],out=os.path.join(work_dir,'_'.join([pID,'SC.png'])))
# for x in pDescs:
# 	print x
# for x in tstPanLst:
# 	print x

#select CPC compounds with the highest average SS and CC
#simple search
work_dir = '/xchip/cogs/projects/ASG_dose_time/plate_design/CPC006/by_SS_CC'
cpcLst = CM.find({'brew_prefix':{'$regex':brewPrfx}},{'sig_id':True,'pert_id':True})
sigLst = [x['sig_id'] for x in cpcLst]
pertLst = [x['pert_id'][:13] for x in cpcLst]
unPerts = set(pertLst)
ssDict = {}
for pID in unPerts:
	CM = mutil.CMapMongo()
	brewPrfx = 'CPC006'
	#query only CPC plates
	cpcLst = CM.find({'brew_prefix':{'$regex':brewPrfx},'pert_id':{'$regex':pID}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
	#query all of db
	# cpcLst = CM.find({'pert_id':{'$regex':pID}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
	ssLst = [x['distil_ss'] for x in cpcLst]
	sigLst = [x['sig_id'] for x in cpcLst]
	cellLst = [x['cell_id'] for x in cpcLst]
	unqCell = set(cellLst)
	nCells = len(unqCell)
	n_sigs = len(sigLst)
	ssDict[pID] = numpy.mean(ssLst)
	if numpy.mean(ssLst) < 6:
		continue
	### make SC plots ###
	# sco = sc.SC()
	# sco.add_sc_from_mongo_sig_ids(sigLst)
	# # sco.set_thresh_by_specificity(0.8)
	# sco.plot(title=(str(n_sigs) + ' signatures of ' + str(pID) + ' in affogato, ' + str(nCells) + 'cell lines'),pos_con=['None'],out=os.path.join(work_dir,'_'.join([pID,'SC.png'])))
	
ssValues = ssDict.values()
ssValues.sort()
#if average ss is greater than about 8
#make ss in each plot


#check for BRD-K11540476
CM = mutil.CMapMongo()
cpcLst = CM.find({'pert_id':{'$regex':'BRD-K11540476'}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
ssLst = [x['distil_ss'] for x in cpcLst]
ccLst = [x['distil_cc_q75'] for x in cpcLst]
sigLst = [x['sig_id'] for x in cpcLst]
cellLst = [x['cell_id'] for x in cpcLst]
unqCell = set(cellLst)
nCells = len(unqCell)
n_sigs = len(sigLst)
sco = sc.SC()
sco.add_sc_from_mongo_sig_ids(sigIDs)
sco.plot(title=('BRD-K11540476'),pos_con=['None'],out=os.path.join(work_dir,'_'.join(['BRD-K11540476','SC.png'])))








### make SC plots of corey's list 
work_dir = '/xchip/cogs/projects/ASG_dose_time/plate_design/potential_list'
fname2 = '/xchip/cogs/projects/ASG_dose_time/plate_design/potential_list/pan_list1.grp'
tstPanLst = open(fname2).read().splitlines()
tstPanLst = [x.split(':')[0] for x in tstPanLst]
pDescs = []
pIDs = []
for pID in tstPanLst:
	CM = mutil.CMapMongo()
	brewPrfx = 'CPC006'
	#query only CPC plates
	cpcLst = CM.find({'brew_prefix':{'$regex':brewPrfx},'pert_id':{'$regex':pID}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
	#query all of db
	# cpcLst = CM.find({'pert_id':{'$regex':pID}},{'sig_id':True,'pert_desc':True,'pert_id':True,'distil_ss':True,'distil_cc_q75':True,'cell_id':True})
	ssLst = [x['distil_ss'] for x in cpcLst]
	sigLst = [x['sig_id'] for x in cpcLst]
	cellLst = [x['cell_id'] for x in cpcLst]
	pDescLst = [x['pert_desc'] for x in cpcLst]
	pDescs.append(pDescLst[0])
	pIDLst = [x['pert_id'] for x in cpcLst]
	pIDs.append(pIDLst[0])
	unqCell = set(cellLst)
	nCells = len(unqCell)
	n_sigs = len(sigLst)
	sco = sc.SC()
	sco.add_sc_from_mongo_sig_ids(sigLst)
	# sco.set_thresh_by_specificity(0.8)
	sco.plot(title=(str(n_sigs) + ' signatures of ' + str(pID) + ' in affogato, ' + str(nCells) + 'cell lines'),pos_con=['None'],out=os.path.join(work_dir,'_'.join([pID,'SC.png'])))
#shorten pIDs
pIDs = [x[:13] for x in pIDs]
pIDset1 = set(pIDs)
#load ASG pert ids
gctfile = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_6H/by_pert_id_pert_dose/ASG001_MCF7_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
db = gct.GCT() #make a gct object
db.read(gctfile)
ASGpertID = db.get_column_meta('pert_id')
ASGpertDesc = db.get_column_meta('pert_desc')
# ASGpertID = [x[:13] for x in ASGpertID]
ASGpIDset = set(ASGpertID)

pIDset1.intersection(ASGpertID)

setDesc = set(ASGpertDesc)
for desc in setDesc:
	i = ASGpertDesc.index(desc)
	pID = ASGpertID[i]
	print desc + ':' + pID

#what are all the trt_cp poscons - make sure to include these

# ### load in ASG plate to get LM genes
# gctfile = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_6H/by_pert_id_pert_dose/ASG001_MCF7_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
# lm = gct.GCT() #make a gct object
# lm.read(gctfile)
# LMrids = lm.get_rids()

# #load in z-score data for signatures from affogato gctx?
# gctfile = '/xchip/cogs/data/build/affogato/affogato_r1_score_n398050x22268.gctx'
# db = gct.GCT() #make a gct object
# db.read(gctfile,cid=sigLst,rid=LMrids)

#mongo test queries:
# dmsoLst = CM.find({'pert_id':{'$regex':qStr},'cell_id':cell1},{'sig_id':True})
# dmsoLst = CM.find({'pert_id':{'$regex':qStr},'cell_id':cell1,'distil_ss':{'$gt':ssMin}},{'sig_id':True})
# dmsoLst = CM.find({'pert_id':{'$regex':qStr},'cell_id':cell1,'distil_ss':{'$gt':ssMin}},{})