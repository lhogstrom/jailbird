#! /usr/bin/env python
'''
frequen itemset analysis:

support = percentage of instances that contain a specific itemset
confidence =  support({ item1, item2})/ support(item1)
--> maybe set(item1, item2) occurs infrequently so has week general support, 
but this set always happens when item1 appears. This is still a meaningful 
finding with high confidence(item1 -> item2) and low support of set(item1,item2)

apriori principal = if an itemset is frequent, then all of its subsets are frequent
or if an itemset is infrequent then its supersets are also infrequen
'''

import cmap.analytics.itemset as itemset
import cmap.io.gct as gct
import cmap.util.mongo_utils as mutil
import numpy as numpy
import pandas as pd
import matplotlib.pyplot as plt

def _all_indices(value, qlist):
	'''
	return all indices where the input string matches the item in the list
	'''
	indices = []
	indx = -1
	while True:
		try:
			indx = qlist.index(value, indx+1)
			indices.append(indx)
		except ValueError:
			break
	# indices = [i for i, x in enumerate(qlist) if x == value]
	return indices

### test itemset 
dataset = itemset._load_sample_dataset()
#since we have a set number of probes, should we modify this somehow?
candSets = itemset.createC1(dataset)
D = map(set,dataset)
L1, support1 = itemset.scanD(D,candSets,.3)
itemset.aprioriGen(L1,2)

### load in dat from gct
work_dir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/scratch'
fname = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
db = gct.GCT()
db.read(fname)


# iMatch = _all_indices('trichostatin-a',db.get_column_meta('pert_desc'))
# selectCids = [db.get_cids()[j] for j in iMatch]
# # generate transactions from the gct object
# up_transactions = itemset.transactions_from_gct(db,cid=selectCids, thresh=2)
# dn_transactions = itemset.transactions_from_gct(db,cid=selectCids, thresh=-2)

# #run the apriori algoritm to generate frequent itemsets
# up_res, supportUp = itemset.apriori(up_transactions, minsupport=.9, max_size=15)
# dn_res, supportDn = itemset.apriori(dn_transactions, minsupport=.5, max_size=15)
# # upSet, dnSet = itemset.gene_sets_from_gct(db,up_thresh=2, dn_thresh=-2, minsupport=0.9, max_size=30)


# ### create a up / dn list based on probe support
LMrid = db.get_rids()
# ### load in data from mongo query
# qStr = 'vorinostat'
# CM = mutil.CMapMongo()
# # cmpdLst = CM.find({'pert_desc':{'$regex':qStr},'cell_id':'PC3'},{'sig_id':True},limit=100)
# cmpdLst = CM.find({'pert_desc':{'$regex':qStr}},{'sig_id':True})
# # cmpdLst = CM.find({'pert_desc':{'$regex':qStr}},{'sig_id':True},limit=20)
# rslt = gct.GCT()
# affName = '/xchip/cogs/data/build/affogato/affogato_r1_score_n398050x22268.gctx'
# rslt.read(affName,cid=cmpdLst,rid=LMrid)
# #find all the probes which pass a specific threshold
# up_transactions = itemset.transactions_from_gct(rslt, thresh=2)
# dn_transactions = itemset.transactions_from_gct(rslt, thresh=-2)
# LenTrans = [len(x) for x in up_transactions]
# #level one analysis
# thresh1 = .5
# C1 = itemset.createC1(up_transactions) #would it be better to just input the 978 probes?
# D = map(set, up_transactions)
# L1, support_data = itemset.scanD(D, C1, thresh1)
# upOut = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/freq_set/up_setLst.txt'
# with open(upOut,'w') as f:
# 	for d in L1:
# 		str1 = ''.join(d) #convert frozen set to string
# 		f.write(str1 + '\n')
# C1 = itemset.createC1(dn_transactions) #would it be better to just input the 978 probes?
# D = map(set, dn_transactions)
# L1, support_data = itemset.scanD(D, C1, thresh1)
# dnOut = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/freq_set/dn_setLst.txt'
# with open(dnOut,'w') as f:
# 	for d in L1:
# 		str1 = ''.join(d) #convert frozen set to string
# 		f.write(str1 + '\n')


#for a series of drugs
# 1) rank order probes acording to number of occurences in training 
# 2) cmap query with a probe list of ~50
# 3) systematically decrease probe list - see how this affects ditribution 
#    of 1) self-connection ES scores 2) query rank distribution


#which compounds have the highest distribution of probe frequency values?
#could do this by: 1) looking at lm_up50 or 2) by doing threshold

qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qDose = db.get_column_meta('pert_dose')

upSupportDict = {}
dnSupportDict = {}
pertSet = set(qPert)
for pert in pertSet:
	iP = _all_indices(pert, qPert) #index of doses on plate
	if len(iP) < 2:
		print pert + ' has only one instance'
		continue
	uDose = [qDose[i] for i in iP]
	fDose = [float(x) for x in uDose] #convert strings to float
	aDose = numpy.asarray(fDose) #convert to numpy array
	iD = aDose.argsort() #local ordering
	sDose = [fDose[j] for j in iD] #sort local doses
	iPo =  [iP[i] for i in iD] #ordered index
	#sMat = ESmat[:,iPo]
	#sMat.sort(axis=0)
	#mongo query for each unique pertID
	qStr = qPertID[iPo[0]] #set pertID
	if len(qStr) >= 13:
		qStr = qStr[0:13] #shorten qPertID
	CM = mutil.CMapMongo()
	#cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True})
	edge50Lst = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True,'up50_lm':True,'dn50_lm':True,'cell_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
	nInstances = len(edge50Lst)
	upProbeDict = {}
	for edge in edge50Lst:
		upPrs = edge['up50_lm']
		for upPr in upPrs:
			if upProbeDict.has_key(upPr):
				count = upProbeDict[upPr][0]
				upProbeDict[upPr][0] = count +1
				upProbeDict[upPr][1] = upProbeDict[upPr][0]/float(len(edge50Lst))
			else:
				upProbeDict[upPr] = [1]
				upProbeDict[upPr].append(1/float(len(edge50Lst)))
	dnProbeDict = {}
	for edge in edge50Lst:
		upPrs = edge['dn50_lm']
		for upPr in upPrs:
			if dnProbeDict.has_key(upPr):
				count = dnProbeDict[upPr][0]
				dnProbeDict[upPr][0] = count +1
				dnProbeDict[upPr][1] = dnProbeDict[upPr][0]/float(len(edge50Lst))
			else:
				dnProbeDict[upPr] = [1]
				dnProbeDict[upPr].append(1/float(len(edge50Lst)))
	upProbeSer = pd.Series(upProbeDict)
	upProbeSer.sort()
	n_top = 50
	#what does the distribution of support scores for n_top look like?
	upSupport = [upProbeSer[-(j+1)][1] for j in range(n_top)]
	upSupportDict[pert] = [len(edge50Lst)] #first write the the number of signatures in the db
	upSupportDict[pert].append(upSupport) #then write the probes with the most support

	#save down support
	dnProbeSer = pd.Series(dnProbeDict)
	dnProbeSer.sort()
	# what does the distribution of support scores for n_top look like?
	dnSupport = [dnProbeSer[-(j+1)][1] for j in range(n_top)]
	dnSupportDict[pert] = [len(edge50Lst)] #first write the the number of signatures in the db
	dnSupportDict[pert].append(dnSupport) #then write the probes with the most support

#graph Support - UP
for pert in upSupportDict:
	if upSupportDict[pert][0] >= 50: #continue only if compound has more than 50 sigs in db
		sprt = upSupportDict[pert][1] #support for top probes for that compound
		plt.plot(sprt)
		plt.ylim(0,1)
plt.xlabel('top probes sorted by support')
plt.ylabel('probe support (percent of db instances)')
plt.title('probe support by compound - up regulated')
plt.show()
#graph Support - DN
for pert in dnSupportDict:
	if dnSupportDict[pert][0] >= 50: #continue only if compound has more than 50 sigs in db
		sprt = dnSupportDict[pert][1] #support for top probes for that compound
		plt.plot(sprt)
		plt.ylim(0,1)
plt.xlabel('top probes sorted by support')
plt.ylabel('probe support (percent of db instances)')
plt.title('probe support by compound - down regulated')
plt.show()

#how does probe support correlate with SS or CC

#L1000_EPSILON.R2.chip is most of affogato -- some instances are based on epsilon-R1 or delta-R1
#does the dict contain non-LM?
notLM = []
for p in upProbeDict:
	if not p in LMrid:
		notLM.append(p)
#does a specific mongo entry contain a non-LM probe id
offSigs = []
for edge in edge50Lst:
	for prb in edge['up50_lm']:
		if not prb in LMrid:
			offSigs.append(edge['sig_id'])

#run the apriori algoritm to generate frequent itemsets
up_res, supportUp = itemset.apriori(up_transactions, minsupport=.9, max_size=10)
dn_res, supportDn = itemset.apriori(dn_transactions, minsupport=.9, max_size=10)



### do comparison different sets as inputs to cmap queries -  which yields higher es scores?
# --> training and test dataset
# how does the length of a query list affect the connection results?
