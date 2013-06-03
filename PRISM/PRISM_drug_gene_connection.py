#! /usr/bin/env python
'''
make a query of dose data
'''

import os
import cmap.io.gct as gct
import glob as glob
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cmap.analytics.dose as doseClass
import subprocess
import time

work_dir = '/xchip/cogs/projects/PRISM/DG_dose_connect'
baseDir = '/xchip/cogs/data/brew/a2y13q1' #where the brewed data lives
outBase = '/xchip/cogs/projects/PRISM/dose_plate_output-by_pert_id_pert_dose'
cellLst = ['A375','MCF7','PC3']
tpLst = ['6H','24H']
plateLst = ['PRISM001']
cntrlT = 'pc'

# ##make system call - run dose_plate_tool
# for cell in cellLst:
# 	outdir = work_dir + '/' + cell
# 	if not os.path.exists(outdir):
# 		os.mkdir(outdir)
# 	#get all CGS for a cell line
# 	CM = mutil.CMapMongo()
# 	CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','cell_id':cell},{'sig_id':True,'pert_iname':True})
# 	if CGSbyCell:
# 		nCGS = len(CGSbyCell)
# 		sigF = os.path.join(outdir, cell+ '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
# 		with open(sigF, 'w') as f:
# 			for sig in CGSbyCell:
# 				f.write(sig['sig_id'] + '\n')
	# for tp in tpLst:
	# 	for platePrefix in plateLst:
	# 		plateName = '_'.join([platePrefix,cell,tp])
	# 		print plateName
	# 		inPut = glob.glob('/'.join([baseDir,
	# 									plateName,
	# 									'by_pert_id_pert_dose',
	# 									plateName+'_COMPZ.MODZ_SCORE_LM_n*x978.gctx']))[0]
# 			dp = doseClass.DosePlate()
# 			dp.add_from_gct(inPut)
# 			dp.examine_doses_tested()			
# 			cellLs = dp.gct.get_column_meta('cell_id')
# 			timePs = dp.gct.get_column_meta('pert_time')
# 			probes = dp.gct.get_rids()
# 			doses = dp.pert_doses
# 			fup = os.path.join(work_dir,cell,cell + '_dose_up_50.gmt')
# 			fdn = os.path.join(work_dir,cell,cell + '_dose_dn_50.gmt')
# 			open(fup,'a') #overwrite existing grp file
# 			open(fdn, 'a') #overwrite existing grp file
# 			n_edge = 50
# 			for pertID in dp.perts_at_dose:
# 				for i in dp.doseIndDict[pertID]: #loop through index of each dose
# 					profile = dp.gct.matrix[:,i]
# 					n_prof = len(profile)
# 					iprofile = profile.argsort() #indices that sort array
# 					iprofile = iprofile[::-1] #switch indicies to decend
# 					sprofile = profile[iprofile]
# 					itop = iprofile[0:(n_edge)]
# 					ibot = iprofile[-n_edge:n_prof]
# 					col_name = pertID + '_' + str(doses[i]) + 'um_' + plateName
# 					ptop = [] 
# 					pbot = []
# 					for j,it in enumerate(itop):
# 						ptop.append(probes[it]) #make probe id list
# 					for j,ip in enumerate(ibot):
# 						pbot.append(probes[ip]) #make probe id list
# 					#write to gmt list 
# 					with open(fup,'a') as f:
# 						f.write(col_name + '\t' + col_name + '\t')
# 						for pt in ptop:
# 							f.write(pt + '\t')
# 						f.write('\n')
# 					with open(fdn,'a') as f:
# 						f.write(col_name + '\t' + col_name + '\t')
# 						for pb in pbot:
# 							f.write(pb + '\t')
# 						f.write('\n')

### query cgs using dose gmts 
# processes = set()
# max_processes = 7
# for cell in cellLst:
# 	cellDir = os.path.join(work_dir,cell) 
# 	cidF = glob.glob(cellDir + '/' + cell + '_all_CGS_sig_ids_n*.grp')
# 	if not cidF:
# 		print 'CGS sig IDs not found'
# 		continue
# 	cidF = cidF[0]
# 	outdir = os.path.join(work_dir,cell,'sig_query_out')
# 	if not os.path.exists(outdir):
# 		os.mkdir(outdir)
# 	fup = os.path.join(work_dir,cell,cell + '_dose_up_50.gmt')
# 	fdn = os.path.join(work_dir,cell,cell + '_dose_dn_50.gmt')
# 	cmd = ' '.join(['rum -q local sig_query_tool',
# 			 '--uptag ' + fup,
# 			 '--dntag ' + fdn,
# 			 '--metric wtcs',
# 			 '--column_space custom',
# 			 '--cid ' + cidF,
# 			 '--out ' + outdir,
# 			 '--mkdir false',
# 			 '--save_tail false'])
# 	# os.system(cmd)
# 	processes.add(subprocess.Popen(cmd,shell=True))
# 	if len(processes) >= max_processes:
# 		os.wait()
# 		processes.difference_update(
# 			p for p in processes if p.poll() is not None)

gp_type = 'KD' # genetic perturbation type
#which cell lines have a result dir
cellDirs = [f for f in os.listdir(work_dir) if os.path.isdir(work_dir+'/'+f)]
prog = progress.DeterminateProgressBar('Drug-target')
df = pd.DataFrame()
dfRank = pd.DataFrame()
#loop through each cell line add to df
# for icell, cell1 in enumerate(cgsCells):
for icell, cell1 in enumerate(cellDirs):
	#define directories and load in outputs
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx'):
		print cell1 + 'no query result file'
		continue #if no results file, skip loop
	rsltFile = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(rsltFile)
	prog.update('analyzing {0}',icell,len(cellDirs))
	rsltF = rslt.frame
	rsltF = rsltF.T
	indVals = rsltF.index.values
	pertVals = [ind[:13] for ind in indVals]
	#make the column name gene and pert time
	geneVals = []
	for ind in rsltF.columns:
		gene = ind.split(':')[1]
		tp = ind.split(':')[0].split('_')[-1]
		gname = '_'.join([gene, tp])
		geneVals.append(gname)
	if len(geneVals) > len(set(geneVals)):
		print 'duplicate CGS for this celline'
	newF = rsltF
	newF.index = [pertVals, rsltF.index.values]
	newF.columns = geneVals
	rankF = newF.rank(ascending=False,axis=1)
	perRankF = rankF / float(rankF.shape[1]) * 100.0
	# newF['cell'] = cell1
	# perRankF['pert'] = pertVals
	# perRankF['cell'] = cell1
	#add cell line result to combined df
	if len(df) == 0:
		df = newF
		dfRank = perRankF
	else:
		df = pd.concat([df,newF],axis=0)
		dfRank = pd.concat([dfRank,perRankF],axis=0)


### get PRISM brds
file1 = '/xchip/cogs/data/brew/a2y13q1/PRISM001_PC3_24H/by_pert_id_pert_dose/PRISM001_PC3_24H_COMPZ.MODZ_SCORE_LM_n60x978.gctx'
pDescDict = {}
dp = doseClass.DosePlate()
dp.add_from_gct(file1)
for i,pert in enumerate(dp.pert_ids):
	if not pDescDict.has_key(pert[:13]):
		pDescDict[pert[:13]] = dp.pert_descs[i]	

### check BEZ235-mtor connection
# inds = dfRank.columns
# # mtors = [ind for ind in inds if ind[:4] == 'MTOR']
# brd = 'BRD-K12184916' #bez235
# # brd = 'BRD-A75409952' #wortmanin
# # BRD-A79768653 #sirolimus 
# # ind = 'MTOR_144H'
# # ind = 'MTOR_96H' #*
# # ind = 'PDK1_96H'
# # ind = 'PIK3CA_144H'
# # ind = 'PIK3CA_96H' #*
# # ind = 'RPTOR_144H'
# ind = 'RPTOR_96H' #*

#gemcitabine
# brd = 'BRD-K15108141' #gemcitabine
# # rrm1s = [ind for ind in inds if ind[:4] == 'RRM1']
# # rrm1s = [ind for ind in inds if ind[:5] == 'ERCC1']
# #check PTEN / ERCC1
# # ind = 'RRM1_144H' #*
# ind = 'RRM1_96H' #*
# ind = 'PTEN_144H'
# ind = 'PTEN_96H' #*
# ind = 'ERCC1_144H'
# ind = 'ERCC1_96H' #*
#check RRM2 and RRM2B

##vorinostat 
# brd = 'BRD-K81418486' #vorinostat 
# # # rrm1s = [ind for ind in inds if ind[:5] == 'HDAC3']
# # # ind = 'HDAC6_144H' 
# # # ind = 'HDAC6_96H'
# # ind = 'HDAC3_144H' #*
# ind = 'HDAC3_96H' 

#teniposide 
# brd = 'BRD-A35588707' #teniposide 
# # rrm1s = [ind for ind in inds if ind[:5] == 'TOP2A']
# ind = 'TOP2A_144H'
# ind = 'TOP2A_96H'

#vx-680
brd = 'BRD-K59369769'
brd = 'BRD-K36740062' #GSK-1070916
brd = 'BRD-K01737880' # cp of interest 
# rrm1s = [ind for ind in inds if ind[:4] == 'AURK']
ind = 'AURKAIP1_144H'
ind = 'AURKAIP1_96H' #maybe
ind = 'AURKA_144H'
ind = 'AURKA_96H' #*
ind = 'AURKB_144H'
ind = 'AURKB_96H'

inFile = '/Users/hogstrom/Desktop/dose_dg_targets.txt'
DrugList = []
cgsList = []
with open(inFile,'rt') as f:
	for string in f:
		splt = string.split('\r')
		for i,line in enumerate(splt):
			if i == 0: # skip headder
				continue
			splt2 = line.split('\t')
			drug1 = splt2[0] #the pert_id listed the line
			cgs = splt2[1]
			DrugList.append(drug1)
			cgsList.append(cgs)

for i,x in enumerate(cgsList):
	brd = DrugList[i]
	ind = cgsList[i]
	cpRank = dfRank.ix[brd]
	cpRes = df.ix[brd]
	rnks = cpRank[ind]
	colName = [x for x in rnks.index]
	dose = [x.split('_')[1][:-2] for x in colName]
	fdose = [float(x) for x in dose]
	logDose = [np.log10(x) for x in fdose]
	#create a list where each dose value increases by 1
	doseSpace = []
	doseSpace.extend(range(1,8))
	doseSpace.extend(range(1,8))
	doseSpace.extend(range(1,8))
	doseSpace.extend(range(1,8))
	doseSpace.extend(range(1,8))
	doseSpace.extend(range(1,8))
	### graph rank
	# plt.scatter(logDose,rnks.values)
	# plt.scatter(fdose,rnks.values)
	# plt.scatter(doseSpace,rnks.values)
	fig1 = plt.figure(figsize = [12,10])
	# ax1 = fig1.add_axes([0.1, 0.1, 0.1, 0.1])
	# h1 = plt.scatter(doseSpace[:14],rnks.values[:14],label='PC3',color='green',s=60,alpha=.4)
	h1 = plt.scatter(doseSpace[14:28],rnks.values[14:28],label='MCF7',color='blue',s=60,alpha=.4)
	h1 = plt.scatter(doseSpace[28:],rnks.values[28:],label='A375',color='purple',s=60,alpha=.4)
	plt.legend()
	plt.ylim((-5,100))
	plt.xlabel('dose um',fontsize=20,fontweight='bold')
	plt.xticks(range(1,8), dose[:7], rotation = 45)
	plt.ylabel('percent rank',fontsize=20,fontweight='bold')
	plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type,fontweight='bold')
	plt.savefig(os.path.join(work_dir,'known_connections',brd +'_' + ind + '_percent_rank.png'))
	plt.close()
	### graph wtcs
	wtcs = cpRes[ind]
	# plt.scatter(logDose,wtcs.values)
	fig1 = plt.figure(figsize = [12,10])
	# h1 = plt.scatter(doseSpace[:14],wtcs.values[:14],label='PC3',color='green',s=60,alpha=.4)
	h1 = plt.scatter(doseSpace[14:28],wtcs.values[14:28],label='MCF7',color='blue',s=60,alpha=.4)
	h1 = plt.scatter(doseSpace[28:],wtcs.values[28:],label='A375',color='purple',s=60,alpha=.4)
	plt.legend()
	plt.ylim((-1,1))
	plt.xlabel('dose um',fontsize=20,fontweight='bold')
	plt.xticks(range(1,8), dose[:7], rotation = 45)
	plt.ylabel('wtcs',fontsize=20,fontweight='bold')
	plt.title(pDescDict[brd] + ' - ' + ind + ' connection - ' + gp_type,fontweight='bold')
	plt.savefig(os.path.join(work_dir,'known_connections',brd +'_' + ind + '_wtcs.png'))
	plt.close()
