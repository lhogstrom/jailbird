#! /usr/bin/env python
#run file like this
#execfile('/xchip/cogs/hogstrom/scripts/dose_timeCourse/ASG_ss_metric.py')

import cmap.analytics.sc as sc
import os.path
import cmap.analytics.signature_strength as ss
import numpy
import scipy
from scipy import stats
import cmap.io.gct as gct
#import ljh_dose_analysis_tool as dose
import cmap.analytics.signature_strength as ss
import cmap.util.mongo_utils as mu

#plot tool
import pylab as pl
import matplotlib.pyplot as plt

#work_dir = '/xchip/cogs/hogstrom/analysis/scratch/Nov27' #MCF7 24h
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/Nov29/dose_analysis_tool.1354211763774' #pc3 6h
cellLine = 'PC3'
timeP = '6H'
gctfile = '/xchip/obelix/pod/brew/pc/ASG001_%s_%s/by_pert_id_pert_dose/ASG001_%s_%s_COMPZ.MODZ_SCORE_LM_n85x978.gctx' % (cellLine,timeP,cellLine,timeP) 

#make a gct object
db = gct.GCT()
db.read(gctfile)

### ss calculations ###
SS1 = ss.SigStrength()
SS1.sig_strength_from_gct_file(gctfile,do_zthresh=False) 
SS2 = ss.SigStrength()
SS2.sig_strength_from_gct_file(gctfile,do_zthresh=True) #ss with threshold

qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qDose = db.get_column_meta('pert_dose')

## plot ss orig with dose
SSin = SS1.ss
ssMax = numpy.nanmax(SSin)
ssMin = numpy.nanmin(SSin)
pearSet = []
spearSet = []
pertSet = set(qPert)
for pert in pertSet:
    iP = all_indices(pert, qPert)
    if len(iP) < 2:
	print pert + ' has only one instance'
	continue
    uDose = [qDose[i] for i in iP]
    uMetric = [SSin[i] for i in iP] #metric drawn from large list
    fDose = [float(x) for x in uDose] #convert strings to float
    aDose = numpy.asarray(fDose) #convert to numpy array
    iD = aDose.argsort()
    sDose = [fDose[j] for j in iD] #sort local doses
    sMetric = [uMetric[j] for j in iD] #sort local metric
    Pear = scipy.stats.pearsonr([1,2,3,4],sMetric)
    Spear = scipy.stats.spearmanr([1,2,3,4],sMetric)
    pearSet.append(Pear)
    spearSet.append(Spear)
    # plot 1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([1,2,3,4],sMetric,'bo-')
    ax.set_ylim(ssMin,ssMax)
    ax.set_xticks([1,2,3,4])
    ax.set_xticklabels( tuple(sDose) )
    ax.set_ylabel('ss')
    ax.set_xlabel('dose')
    ax.set_title(pert + ' signature stength')
    fname = pert + '_ss.png'
    outf = os.path.join(work_dir,fname)
    plt.savefig(outf, bbox_inches=0)

#plot ss with thresh with dose##
SSin = SS2.ss
ssMax = numpy.nanmax(SSin)
ssMin = numpy.nanmin(SSin)
pearSet = []
spearSet = []
pertSet = set(qPert)
for pert in pertSet:
    iP = all_indices(pert, qPert)
    if len(iP) < 2:
	print pert + ' has only one instance'
	continue
    uDose = [qDose[i] for i in iP]
    uMetric = [SSin[i] for i in iP] #metric drawn from large list
    fDose = [float(x) for x in uDose] #convert strings to float
    aDose = numpy.asarray(fDose) #convert to numpy array
    iD = aDose.argsort()
    sDose = [fDose[j] for j in iD] #sort local doses
    sMetric = [uMetric[j] for j in iD] #sort local metric
    Pear = scipy.stats.pearsonr([1,2,3,4],sMetric)
    Spear = scipy.stats.spearmanr([1,2,3,4],sMetric)
    pearSet.append(Pear)
    spearSet.append(Spear)
    # plot 1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([1,2,3,4],sMetric,'ro-')
    ax.set_ylim(ssMin,ssMax)
    ax.set_xticks([1,2,3,4])
    ax.set_xticklabels( tuple(sDose) )
    ax.set_ylabel('ss')
    ax.set_xlabel('dose')
    ax.set_title(pert + ' signature stength Z THRESH')
    fname = pert + '_ss_thresh.png'
    outf = os.path.join(work_dir,fname)
    plt.savefig(outf, bbox_inches=0)

## make qq plot
ESmat = db.matrix
col1 = ESmat[:,0]
#calculate null distribution
mu, sigma = 0, 1
s = numpy.random.normal(mu, sigma, len(ESmat[:,1]))
s.sort()

#qq plot
pertSet = set(qPert)
for pert in pertSet:
    iP = _all_indices(pert, qPert) #index of doses on plate
    if len(iP) < 2:
	print pert + ' has only one instance'
	continue
    uDose = [qDose[i] for i in iP]
    uMat = ESmat[:,iP]
    fDose = [float(x) for x in uDose] #convert strings to float
    aDose = numpy.asarray(fDose) #convert to numpy array
    iD = aDose.argsort()
    sDose = [fDose[j] for j in iD] #sort local doses
    sMat = uMat[:,iD] #sort local metric
    sMat.sort(axis=0)
    #make qq plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(s,s,'b')
    ax.plot(s,sMat[:,0],'r.',label="dose = .08um")
    ax.plot(s,sMat[:,1],'g.',label="dose = .4um")
    ax.plot(s,sMat[:,2],'y.',label="dose = 2um")
    ax.plot(s,sMat[:,3],'c.',label="dose = 10um")
    ax.legend(loc=4)
    #plt.show()
    ax.set_ylabel('observed z-score')
    ax.set_xlabel('expected z-score')
    ax.set_title(pert + ' signature stength Z THRESH')
    fname = pert + '_qq.png'
    outf = os.path.join(work_dir,fname)
    plt.savefig(outf, bbox_inches=0)


#ax.plot(s,sMat[:,0],'.',markersize=3)

#qPert = db.get_column_meta('pert_desc')
#qPertID = db.get_column_meta('pert_id')
#qDose = db.get_column_meta('pert_dose')
#probeIDs = db.get_row_meta('id')

### retrieve enrichment info for each probe
#qStr = qPertID[i]
#cmpd1 = qPert[i]
#dose1 = qDose[i]
#if len(qStr) >= 13:
	#qStr = qStr[0:13] #shorten qPertID

#CM = mu.CMapMongo()
##cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True})
#edge50Lst = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True,'up50_lm':True,'dn50_lm':True,'cell_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
#probeCnts = [0] * len(probeIDs)
#for j,inst in enumerate(edge50Lst):
	#up50 = edge50Lst[j]['up50_lm']
	##loop through every gene in the top and bottom list - where does it live on the rank list?
	#for prb in up50:
		#if prb in probeIDs:
			#iPrb = probeIDs.index(prb)
			#probeCnts[iPrb] = probeCnts[iPrb] +1
	##loop through each dose
#sprobeCnts = [probeCnts[l] for l in iD] #sort probe counts acording to z-score


#g1 = set(probeIDs)
#g2 = set(up50)
#g2.difference(g1)

#Out[359]: 
#set([u'202032_s_at',
     #u'202798_at',
     #u'203562_at',
     #u'204925_at',
     #u'205548_s_at',
     #u'209911_x_at',
     #u'212726_at',
     #u'218255_s_at',
     #u'221848_at'])
#Are these LM probes?


#the badass qq plot
qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qDose = db.get_column_meta('pert_dose')
probeIDs = db.get_row_meta('id')
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
	CM = mu.CMapMongo()
	#cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True})
	edge50Lst = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True,'up50_lm':True,'dn50_lm':True,'cell_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
	#count number of times a probe is in the top/bottom 50 genes of an instance
	upProbeCnts = [0] * len(probeIDs)
	dnProbeCnts = [0] * len(probeIDs)
	for j,inst in enumerate(edge50Lst):
		up50 = edge50Lst[j]['up50_lm']
		dn50 = edge50Lst[j]['dn50_lm']
		#loop through every gene in the top and bottom list - where does it live on the rank list?
		for prb in up50:
			if prb in probeIDs:
				iPrb = probeIDs.index(prb)
				upProbeCnts[iPrb] = upProbeCnts[iPrb] +1
		for prb in dn50:
			if prb in probeIDs:
				iPrb = probeIDs.index(prb)
				dnProbeCnts[iPrb] = dnProbeCnts[iPrb] +1
	#loop through each dose
	for d in iPo:
	#count probe enrichment and plot
			cmpd1 = qPert[d]
			dose1 = qDose[d]
			zLst = ESmat[:,d]
			iLst = zLst.argsort() #sort z-scores and save index
			sLst = zLst[iLst]
			sUpProbeCnts = [upProbeCnts[l] for l in iLst] #sort probe counts acording to z-score
			sDnProbeCnts = [dnProbeCnts[l] for l in iLst]
			#mkrs = numpy.sqrt(sprobeCnts) # non linear scaling of marker points
			sUpProbeCnts = [float(l) for l in sUpProbeCnts] #convert to float
			sDnProbeCnts = [float(l) for l in sDnProbeCnts] #convert to float
			upPercMkrs = numpy.divide(sUpProbeCnts,max(sUpProbeCnts)) #divide by max count to make for relative frequency
			dnPercMkrs = numpy.divide(sDnProbeCnts,max(sDnProbeCnts))
			upMkrs = numpy.multiply(upPercMkrs,100)
			dnMkrs = numpy.multiply(dnPercMkrs,100)
			fig = plt.figure()
			ax = fig.add_subplot(111)
			ax.plot(s,s,'b')
			for j,sl in enumerate(sLst):
				ax.plot(s[j],sl,'r.',markersize=upMkrs[j],alpha=.25)
				ax.plot(s[j],sl,'b.',markersize=dnMkrs[j],alpha=.25)
			ax.set_ylabel('observed z-score')
			ax.set_xlabel('expected z-score')
			ax.set_title(pert + ' signature stength Z THRESH, dose = ' + dose1)
			fname = pert + '_' + dose1 + 'um_internal-external_qq.png'
			outf = os.path.join(work_dir,fname)
			plt.savefig(outf, bbox_inches=0)
