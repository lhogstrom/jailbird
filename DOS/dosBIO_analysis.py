#! /usr/bin/env python
'''
analyze the DOSBIO plates - make sc plots etc
'''

import os
import cmap.io.gct as gct
import cmap.analytics.sc as sc
import glob as glob
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plate = 'DOSBIO001'
cellLine = 'PC3'
refControl ='pc'
timeP = '24H'
dataDir = '/xchip/cogs/projects/DOS/data'
gctfile = glob.glob(dataDir +'/%s/%s_%s_%s/by_pert_id_pert_dose/%s_%s_%s_COMPZ.MODZ_SCORE_LM_*x978.gctx' % (refControl,plate,cellLine,timeP,plate,cellLine,timeP))
gctfile = gctfile[0]

work_dir = '/xchip/cogs/projects/DOS/DOSBIO'
# work_dir = '/xchip/cogs/hogstrom/analysis/scratch/prism/%s_%s_%s/fulcrum' % (cell,timeP,refControl)
# if not os.path.exists(work_dir):
# 	os.mkdir(work_dir)

gcto=gct.GCT()
gcto.read(gctfile)

sco = sc.SC()
sco.add_sc_from_gctx_meta(gctfile, verbose=False)
sco.set_thresh_by_specificity(0.8)
t1 = '%s_%s_%s' % (plate,cellLine,timeP)
outF = os.path.join(work_dir,'_'.join([t1,'SC.png']))
sco.plot(title=t1,out=outF)


pDescs = gcto.get_column_meta('pert_desc')
iBRDs = [i for i, x in enumerate(pDescs) if x[:3] == 'BRD']
ss= gcto.get_column_meta('distil_ss')
ssBRD = [ss[i] for i in iBRDs]
maxSS = max(ssBRD) 
iMaxSS = [i for i,x in enumerate(ss) if x == maxSS][0]
maxBRD = pDescs[iMaxSS]

brdCounts = []
fullPertList = []
for ipert in iBRDs:
	pert = pDescs[ipert]
	CM = mutil.CMapMongo()
	pert_List = CM.find({'pert_id':{'$regex':pert}},{'sig_id':True,'cell_id':True})
	if pert_List:
		brdCounts.append(len(pert_List))
		fullPertList.extend(pert_List)

cellsAll = [sig['cell_id'] for sig in fullPertList]
uniqCells = list(set(cellsAll))
### make gmt signature of drugs of interest
for cell1 in uniqCells:
	outdir = os.path.join(work_dir,cell1)
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	sigIDlist = []
	for sig in fullPertList:
		if sig['cell_id'] == cell1:
			sigIDlist.append(sig['sig_id'])
	sigIDlist = list(set(sigIDlist))
	#write drug signatures by cell line to a file
	sigF = os.path.join(outdir,cell1 + '_cp_sig_ids.grp')
	with open(sigF, 'w') as f:
		[f.write(x + '\n') for x in sigIDlist]
	#get all CGS for a cell line
	CM = mutil.CMapMongo()
	CGSbyCell = CM.find({'pert_type':'trt_sh.cgs','cell_id':cell1},{'sig_id':True,'pert_iname':True})
	if CGSbyCell:
		nCGS = len(CGSbyCell)
		sigF = os.path.join(outdir, cell1+ '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
		with open(sigF, 'w') as f:
			for sig in CGSbyCell:
				f.write(sig['sig_id'] + '\n')

for cell1 in uniqCells:
	cellDir = os.path.join(work_dir,cell1) 
	cidF = glob.glob(cellDir + '/' + cell1 + '_all_CGS_sig_ids_n*.grp')
	if not cidF:
		continue
	cidF = cidF[0]
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not os.path.exists(outdir):
		os.mkdir(outdir)
	# sigF = os.path.join(cellDir, cell1 + '_all_CGS_sig_ids_n' + str(nCGS) + '.grp')
	sigF = os.path.join(cellDir,cell1 + '_cp_sig_ids.grp')
	cmd = ' '.join(['rum -q local sig_query_tool',
			 '--sig_id ' + sigF,
			 '--metric wtcs',
			 '--column_space custom',
			 '--cid ' + cidF,
			 '--out ' + outdir,
			 '--mkdir false',
			 '--save_tail false'])
	os.system(cmd)

prog = progress.DeterminateProgressBar('Drug-target')
for icell, cell1 in enumerate(uniqCells):
	celldir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx'):
		print cell1 + 'no query result file'
		continue #if no results file, skip loop
	rsltF = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(rsltF)
	prog.update('analyzing {0}',icell,len(uniqCells))
	queryRids = rslt.get_rids()
	# for KD:
	queryGenes = [x.split(':')[1] for x in queryRids]
	# for OE:
	# queryGenes = [x.split('_')[1] for x in queryRids]
	resCids = rslt.get_cids()
	resPerts = [x.split(':')[1] for x in resCids]
	#build rank matrix
	sortMatrix = np.argsort(rslt.matrix,axis=0)[::-1]
	rankMatrix = np.zeros_like(sortMatrix)
	for icol in range(rankMatrix.shape[1]):
		col = sortMatrix[:,icol]
		rankMatrix[col,icol] = np.arange(len(sortMatrix[:,icol])) + 1
	#pandas rank
	pRslt = rslt.frame
	rankFrame = pRslt.rank(ascending=False)

	for gene in list(set(queryGenes)):
		for pert in list(set(resPerts)):
			iperts = [i for i,x in enumerate(resPerts) if x[:13] == pert]
			for instance, ipert in enumerate(iperts):
				instance = 

rFrame = rankFrame

geneList = []
#build hierarchical index
for rid in rankFrame.index:
	gene = rid.split(':')[1]
	geneList.append(gene)



rFrame.ix['CGS001_A375_96H:A2M:1','DOS054_A375_24H:BRD-K42543764-001-01-8:5']
hFrame.index.names = ['gene','cellLine','instance']

### strategy 1 - every pairwise comparison is its own row
df = pd.DataFrame()
prog = progress.DeterminateProgressBar('Drug-target')
for icell, cell1 in enumerate(uniqCells):
	celldir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx'):
		print cell1 + 'no query result file'
		continue #if no results file, skip loop
	rsltF = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(rsltF)
	prog.update('analyzing {0}',icell,len(uniqCells))
	flatSeries = rslt.frame.unstack()
	flatFrame = pd.DataFrame(flatSeries,columns=['wtcs'])
	flatFrame.index.names = ['sig_id', 'cgs']
	flatFrame['cell'] = cell1
	indVals = flatFrame.index.values
	pertVals = [ind[0].split(':')[1][:13] for ind in indVals]
	geneVals = [ind[1].split(':')[1] for ind in indVals]
	flatFrame['pert'] = pertVals
	flatFrame['gene'] = geneVals
	df.append(flatFrame)

indLim1 = df['pert'] == 'BRD-K70875408'
indLim2 = df['gene'] == 'BBS9'
df[indLim1 & indLim2]




### strategy 2:
# Create a pandas dataframe that lets you see connection results across 
# cell lines it is structured as follows:
# 	index1 = BRD short
# 	index2 = perurbation sig_id
# 	each column - a unique gene ID/ time point - representing the CGS for that gene, matching cell line
# 	cell line listed as a column
work_dir = '/xchip/cogs/projects/DOS/DOSBIO'
#which cell lines have a result dir
cellDirs = [f for f in os.listdir(work_dir) if os.path.isdir(work_dir+'/'+f)]
prog = progress.DeterminateProgressBar('Drug-target')
df = pd.DataFrame()
#loop through each cell line add to df
for icell, cell1 in enumerate(cellDirs):
	#define directories and load in outputs
	celldir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx'):
		print cell1 + 'no query result file'
		continue #if no results file, skip loop
	rsltF = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(rsltF)
	prog.update('analyzing {0}',icell,len(cellDirs))
	rsltF = rslt.frame
	rsltF = rsltF.T
	indVals = rsltF.index.values
	pertVals = [ind.split(':')[1][:13] for ind in indVals]
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
	newF['pert'] = pertVals
	newF['cell'] = cell1
	#add cell line result to combined df
	if len(df) == 0:
		df = newF
	else:
		df = pd.concat([df,newF],axis=0)

#get unique brds with connection data
fullBRDs = []
for ind in df.index:
	brd = ind[0]
	fullBRDs.append(brd)
uniqBRDs = list(set(fullBRDs))
nCPs = []
# for brd in uniqBRDs:
iMany = [i for i,x in enumerate(nCPs) if x > 12]
for ipert in iMany:
	brd = uniqBRDs[ipert]
	cpRes = df.ix[brd]
	nCPs.append(cpRes.shape[0])
	# cntList = []
	# find the number of drug-target instances
	# for col in df.columns:
	# 	if col == 'pert' or col == 'cell':
	# 		print col
	# 		continue
	# 	else:
	# 		mtch = cpRes[col]
	# 		mtchSum = mtch.sum()
	# 		cnts = len(pd.value_counts(mtch))
	# 		if type(cnts) == int:
	# 			cntList.append(cnts)
	meanSer = cpRes.mean()
	# meanVals = meanSer.values
	# ordMean = meanSer.order().index
	# cpRes[ordMean]
	nullCnt = pd.isnull(cpRes)
	#how many cell lines were both the pert and target tested in
	valCounts = nullCnt.shape[0] - nullCnt.sum(axis=0)
	CS_thresh = .45 #theshold for mean ss
	for ind in meanSer[meanSer > CS_thresh].index:
		# print ind
		if valCounts[ind] > 4:
			print brd + ' ' + ind
			#cs wadden gram
			cpRes[ind]
			sKeysStr = []
			for i,cs in enumerate(cpRes[ind]):
				if pd.isnull(cs):
					continue
				else:
					sKeysStr.append(cpRes.index[i])
					# percRanks = targetRnkPercs[cell1][pert][target]
					# yVals = np.repeat(i+1,len(cs))
					yVals = i+1
					# plt.scatter(percRanks,yVals)
					plt.scatter(cs,yVals)
			plt.xlim((-1, 1))
			plt.ylim((0,i+2))
			plt.yticks(range(1, i + 2), sKeysStr, rotation = 0)
			plt.xlabel('wtcs')
			plt.ylabel('pert')
			plt.title(brd + ' - ' + ind + ' connection')
			plt.savefig(os.path.join(work_dir,brd +'_' + ind + '_connections.png'))
			plt.close()

for cs in cpRes[ind]:
	if pd.isnull(cs):
		continue
	else:
		print cs

#get brds where we tested in a bunch of cell lines:




	iSortmean = meanSer.argsort()
	iSort = iSortmean.values
	iSort = iSort.astype(int)
	meanSer[iSort]


#mimic df on a smaller scale:
data1 = pd.Series(np.random.randn(10),index=[['a', 'a', 'a', 'b', 'b', 'b', 'c', 'c', 'd', 'd'],[1, 2, 3, 1, 2, 3, 1, 2, 2, 3]])
data2 = pd.Series(np.random.randn(4),index=[['added1', 'jump', 'cap', 'blot'],[1, 2, 3, 1]])
dd = pd.DataFrame(np.random.randn(10, 4),columns=['g','b','y','q'])
d1 = dd[:3]
d1.index=['a1','a2','a3']
d1 = d1[['b','y']]
d2 = dd[3:7]
d2.index = ['b33', 'b34', 'b35','b36']
pd.merge(d1,d2,how='outer')
pd.concat([d1,d2],axis=0)
pieces = [dd[:3], dd[3:7], dd[7:]]

#make small frame
#be sure to deal with duplicate columns
g = df.ix[:3]
g= g[['A2M_96H','ZRSR2_144H']]
f = df.ix[5:7]
# f= f[['A2M_96H','ZRSR2_144H']]
f= f[['A2M_96H','ZRSR2_144H','A2M_96H']]
# f= f[['ABCB1','ABCB4']]
pd.concat([f,g],axis=0,join='outer')
pd.concat([f,g],axis=0,join_axes)
pd.merge(g,f,how='outer')

#re-index
g.index = ['BRD-K12683703','BRD-K03050720','BRD-K16625384']
f.index = ['BRD-K79734497','BRD-K57200229']

#find drug duplicates
[x for i, x in enumerate(pertVals) if pertVals.count(x) > 1]
[x for i, x in enumerate(geneVals) if geneVals.count(x) > 1]
#get index of one CGS dup
idup = [i for i, x in enumerate(geneVals) if x == 'ZNF134']
rslt.frame.index[idup]
#re-index a series
k4Series = rFrame['DOS054_A375_24H:BRD-K42543764-001-01-8:5']
newInd = []
for ind in k4Series.index:
	newInd.append(ind.split(':')[1])
knew = k4Series.reindex(newInd)
knew = pd.Series(k4Series.values)




ind1 = []
ind2 = []
for cell1 in uniqCells:
	for num in range(4):
		ind1.append(cell1)
		ind2.append(num)

hSeriess = pd.Series(np.zeros_like(ind2),index=[ind1,ind2])
hSeriess = pd.Series(range(len(uniqCells)),index=[uniqCells])

