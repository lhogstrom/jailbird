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
import cmap.analytics.signature_strength as ss
import cmap.util.mongo_utils as mu
# from collections import Counter

#plot tool
import pylab as pl
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
    return indices


#work_dir = '/xchip/cogs/hogstrom/analysis/scratch/Nov27' #MCF7 24h
# work_dir = '/xchip/cogs/hogstrom/analysis/scratch/Nov29/dose_analysis_tool.1354211763774' #pc3 6h
work_dir = '/xchip/cogs/hogstrom/analysis/scratch/ASG_qqPlots/test1'
cellLine = 'MCF7'
timeP = '24H'
gctfile = '/xchip/obelix/pod/brew/pc/ASG001_%s_%s/by_pert_id_pert_dose/ASG001_%s_%s_COMPZ.MODZ_SCORE_LM_n85x978.gctx' % (cellLine,timeP,cellLine,timeP) 

#make a gct object
db = gct.GCT()
db.read(gctfile)

## load in results
frslt = '/xchip/cogs/hogstrom/analysis/scratch/Nov20/dose_analysis_tool.1353449771597/nov20/my_analysis.query_tool.2012112017162991/result_ESLM.COMBINED_n85x398050.gctx'
rslt = gct.GCT()
rslt.read(frslt)
rSigIds = rslt.get_rids()

## make qq plot
ESmat = db.matrix
col1 = ESmat[:,0]
#calculate null distribution
mu1, sigma = 0, 1
s = numpy.random.normal(mu1, sigma, len(ESmat[:,1]))
s.sort()

#the badass qq plot
qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qDose = db.get_column_meta('pert_dose')
probeIDs = db.get_row_meta('id')
# pertSet = set(qPert)
pertSet = [qPert[84]]
for pert in pertSet:
    # pert = qPert[0]
    iP = _all_indices(pert, qPert) #index of doses on plate
    if len(iP) < 2:
    	print pert + ' has only one instance'
    	continue
    uDose = [qDose[i] for i in iP]
    fDose = [float(x) for x in uDose] #convert strings to float
    aDose = numpy.asarray(fDose) #convert to numpy array
    iD = aDose.argsort() #local ordering
    sDose = [fDose[j] for j in iD] #sort local doses
    iPo =  [iP[i] for i in iD] #ordered index acording to dose
    d = 3 #dose of interest - highest?
    rowQ = rslt.matrix[:,iPo[d]] #sort affogato results 
    #mongo query for each unique pertID
    irowQ = rowQ.argsort()[::-1] #sort results 
    SigQ = [rSigIds[irowQ[0]]] #sig ID - top/bottom of ranking profile 
	# qStr = qPertID[iPo[0]] #set pertID as SELF CONNECTION
    anntQ = mu.sig_info(SigQ)
    pertQ = anntQ['pert_desc']
    qStr = anntQ['pert_id'] #
    if len(qStr) >= 13:
        Str = qStr[0:13] #shorten qPertID
    if pertQ == -666.0: #BRD compounds have no pert_desc - fill with pertID
        pertQ = qStr
    CM = mu.CMapMongo()
    #cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True})
    edge50Lst = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True,'up50_lm':True,'dn50_lm':True,'cell_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
    nInstances = len(edge50Lst) #number of instances in db
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
            # upPercMkrs = numpy.divide(sUpProbeCnts,max(sUpProbeCnts)) #divide by max count to make for relative frequency
            # dnPercMkrs = numpy.divide(sDnProbeCnts,max(sDnProbeCnts))
            upPercMkrs = numpy.divide(sUpProbeCnts,nInstances) #divide by total instances to make for relative frequency
            dnPercMkrs = numpy.divide(sDnProbeCnts,nInstances)            
            upMkrs = numpy.multiply(upPercMkrs,1000)
            dnMkrs = numpy.multiply(dnPercMkrs,1000)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(s,s,'b')
            g1 = ax.scatter(s,sLst,upMkrs,c='r',alpha=.10)
            g2 = ax.scatter(s,sLst,dnMkrs,c='b',alpha=.10)
            ax.set_ylabel('observed z-score')
            ax.set_xlabel('expected z-score')
            #set legend based on the number of 
            # r1 = ax.plot(0,0,'r.',markersize=100,alpha=.25)
            # b1 = ax.plot(0,0,'b.',markersize=100,alpha=.25)
            legStrUp = 'probe in 100% of ' + str(nInstances) + ' UP instances'
            legStrDn = 'probe in 100% of ' + str(nInstances) + ' DN instances'
            plt.legend([g1, g2], [legStrUp, legStrDn], numpoints=1, loc=4)
            # plt.legdend([b1], ['probe in 100% of ' + str(nInstances) + 'instances' ], numpoints=1)
            ax.set_title(pert + ' dose = ' + dose1 + 'connection w/ ' + pertQ)
            fname = pert + '_' + dose1 + 'um_' + pertQ + 'um_connection_qq.png'
            outf = os.path.join(work_dir,fname)
            plt.savefig(outf, bbox_inches=0)


#pick the top 100 compounds from the result gctx - which ones are unique?
nt = 100 #number of ranks
SigQs = [rSigIds[irowQ[y]] for y in range(nt)]
anntQs = mu.sig_info(SigQs)
pertQs = [anntQs[y]['pert_desc'] for y in range(nt)]
pertIDQs = [anntQs[y]['pert_id'][0:13] for y in range(nt)] #shorten to 13 characters
#identify which pertIDs occured most often
IDcnts = _get_counts(pertIDQs)
#sort by to occuring pertIDs
value_key_pairs = [(count, Id) for Id, count in IDcnts.items()]
value_key_pairs.sort(reverse=True)

#then pick out all their instances in affogato - list pertID:sig1,sig2,sig4





## affogato perturbation list
fname1 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_desc__affogato.txt';
affPert = open(fname1).read().splitlines()
fname2 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_id_affogato.txt';
affPertID = open(fname2).read().splitlines()
fname3 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/sig_ids_affogato.txt';
affSigID = open(fname3).read().splitlines()
#short pertIDs for compounds
srtLst = []
for ID in affPertID:
    if ID[0:3] == 'BRD' and len(ID) >= 13:
        ID = ID[0:13] #shorten qPertID
        srtLst.append(ID)
    else:
        srtLst.append(ID)
#connect pert annotations
annts = zip (srtLst,affPert,affSigID)
anntCnts = _get_counts(annts,0) #count the instances of each pertID
value_key_pairs = [(count, Id) for Id, count in IDcnts.items()]
value_key_pairs.sort(reverse=True)
pertDict = _get_counts2(annts) #group the sigIDs by pertIDs

#for a unique pertID, find instances in rank file
ESmat = rslt.matrix
iES = ESmat.argsort(axis=0)[::-1] #sort ascending
i = 4 #which column in the matrix
iE = iES[:,i] #ES sort index for one column
sSigID = []
for y in iE:
    sSigID.append(rsltSigID[y]) #make sorted sig ID list

cmpd1 = anntCnts[]
cmpdSigIds = pertDict[cmpd1]    
i1 = [sSigID.index(y) for y in cmpdSigIds] #where instances of the compound of interest sit on the rank list
if len(i1) < 1:
    print cmpd1 + ' has no instances in the cmap database'
    continue
i2 = numpy.array(i1) #convert list to numpy array


def _get_counts(sequence,indx):
    counts = {}
    for x in sequence:
        if x[indx] in counts:
            counts[x[indx]] += 1
        else:
            counts[x[indx]] = 1
    return counts

def _get_counts2(sequence):
    counts = {}
    for ID,pert,sigID in sequence:
        counts.setdefault(ID, []).append(sigID)
    return counts
