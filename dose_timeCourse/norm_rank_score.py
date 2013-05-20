#! /usr/bin/env python
#run file like this
#execfile('/xchip/cogs/hogstrom/scripts/dose_timeCourse/norm_rank_score.py')

import cmap.analytics.sc as sc
import os.path
import cmap.analytics.signature_strength as ss
import numpy
import scipy
from scipy import stats
import cmap.io.gct as gct
#import ljh_dose_analysis_tool as dose

#plot tool
import pylab as pl
import matplotlib.pyplot as plt

work_dir = '/xchip/cogs/hogstrom/analysis/scratch/Nov29/dose_analysis_tool.1354211763774'
cellLine = 'PC3'
timeP = '6H'
gctfile = '/xchip/obelix/pod/brew/pc/ASG001_%s_%s/by_pert_id_pert_dose/ASG001_%s_%s_COMPZ.MODZ_SCORE_LM_n85x978.gctx' % (cellLine,timeP,cellLine,timeP) 

#make a gct object
db = gct.GCT()
db.read(gctfile)

## load in results
#frslt = '/xchip/cogs/hogstrom/analysis/scratch/Nov20/dose_analysis_tool.1353449771597/nov20/my_analysis.query_tool.2012112017162991/result_ESLM.COMBINED_n85x398050.gctx' #mcf7 24h
#frslt = '/xchip/cogs/hogstrom/analysis/scratch/Nov29/dose_analysis_tool.1354211494574/nov29/my_analysis.query_tool.2012112912515791/result_ESLM.COMBINED_n85x398050.gctx' #mcf7 6h
frslt = '/xchip/cogs/hogstrom/analysis/scratch/Nov29/dose_analysis_tool.1354211763774/nov29/my_analysis.query_tool.2012112912562991/result_ESLM.COMBINED_n85x398050.gctx' #pc3 6h
rslt = gct.GCT()
rslt.read(frslt)

fname1 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_desc__affogato.txt';
affPert = open(fname1).read().splitlines()
fname2 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_id_affogato.txt';
affPertID = open(fname2).read().splitlines()
fname3 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/sig_ids_affogato.txt';
affSigID = open(fname3).read().splitlines()
rsltSigID = rslt.get_rids()

#match pertID and pertDesc with sigID for results
rsltPert = []
rsltPertID = []
for i, sig in enumerate(affSigID):
    i1 = affSigID.index(sig)
    rsltPert.append(affPert[i1])
    rsltPertID.append(affPertID[i1])

#write lists to a file
fname1 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_desc_affRslt.txt';
with open(fname1,'a') as f:
    for r in rsltPert:
        f.write(r + '\n')
fname2 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_id_affRslt.txt';
with open(fname1,'a') as f:
    for r in rsltPertID:
        f.write(r + '\n')
fname3 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/sig_ids_affRslt.txt';
with open(fname1,'a') as f:
    for r in rsltSigID:
        f.write(r + '\n')


### read in pertIDs and pertDescs that match the query result gctx ###
fname1 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_desc_affRslt.txt';
rsltPert = open(fname1).read().splitlines()
fname2 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_id_affRslt.txt';
rsltPertID = open(fname2).read().splitlines()
rsltSigID = rslt.get_rids()



#create dictionary of pertIDs
idDict = {} # Here the value is a list: first item is the pert_desc of the key, the second item is a list of indices for that perturbation
for i,x in enumerate(rsltPertID):
    if x[:4] == 'BRD-':
        x = x[:13]
    if idDict.has_key(x):
        idDict[x][1].append(i)
    else:
        pDesc = rsltPert[i]
        idDict[x] = [pDesc]
        idDict[x].append([i])

#check to see if indices are calling correctly 
test1 = 'trichostatin A'
# test1 = 'wortmannin'
i1 = rsltPert.index(test1)
pID1 = rsltPertID[i1]
if pID1[:4] == 'BRD-':
    pID1 = pID1[:13]
inds1 = idDict[pID1][1]


# sort data matrix
ESmat = rslt.matrix #order matrix acording to sigID sort
iES = ESmat.argsort(axis=0)[::-1] #sort ascending 
n_inst = len(iES[:,1])

#find self connections events using dictionary
qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qDose = db.get_column_meta('pert_dose')
qID = db.get_cids()
rPertID = rslt.get_cids()

for i, x in enumerate(qPert):
    iE = iES[:,i] #ES sort index for one column
    qPertDesc = qPert[i]
    qStr = qPertID[i] #was the perturbation being queried
    if qStr[:4] == 'BRD-':
        qStr = qStr[:13] #shorten long BRD ids
    rPertDesc = idDict[qStr][0]
    if qPertDesc != rPertDesc:
        print 'query and result pert_desc do not match'
    iInstances = idDict[qStr][1]
    rnkInsts = [iE[i] for i in iInstances]



#sort sigIDs between query result and affogato saved
iaffSigID = sorted(range(len(affSigID)), key = affSigID.__getitem__) #index of sorting
irsltSigID = sorted(range(len(rsltSigID)), key = rsltSigID.__getitem__) #index of sorting

#list to be generated according to sigID sort index
sRid = []
sAid = []
sApert = []
sApertID = []
for x in irsltSigID:
    sRid.append(rsltSigID[x]) #make sorted sig ID list

for x in iaffSigID:
    sAid.append(affSigID[x]) #make sorted sig ID list
    sApert.append(affPert[x]) #sort perts 
    sApertID.append(affPertID[x]) #sort pertIDs 

if not sAid == sRid:
    print 'result and affogato sig ids are not the same when sorted'



#sortESmat = ESmat[iES]

#loop through each of the perts
avRnk = []
medRnk = []
for i, x in enumerate(qPert):
    iE = iES[:,i] #ES sort index for one column
    IDsorted = []
    for j in iE:
		IDshort = sApertID[j]
		IDshort = IDshort[0:13]
		IDsorted.append(IDshort) 
    qStr = qPertID[i]
    cmpd1 = x
    dose1 = qDose[i]
    if len(qStr) >= 13:
        qStr = qStr[0:13] #shorten qPertID
    #i1 = IDsorted.index(qStr) #give first index of match
    i1 = all_indices(qStr,IDsorted)
    i2 = numpy.array(i1) #convert list to numpy array
    avr = sum(i2)/len(i2) #what is the average ES rank
    md = numpy.median(i2) # what is the median ES rank
    nAv = float(avr)/n_inst #normalize acording to number of instances in db
    nMd = float(md)/len(iES[:,1]) #normalized median
    avRnk.append(nAv) #store average ES rank
    medRnk.append(nMd)
    #plot 
    fname = cmpd1 + '_' + dose1 + '_query_rank.png'
    outf = os.path.join(work_dir,fname)
    fig = plt.figure(figsize=(8.0, 2.0))
    ax = fig.add_subplot(111)
    # the histogram of the data
    n, bins, patches = ax.hist(i2, 30, facecolor='green', alpha=0.75)
    #ax.set_xlim(0, n_inst)
    ax.set_xlim(0, int(round(n_inst,-5))) #round instances to nearest 100k
    ax.set_xlabel('query rank')
    ax.set_ylabel('freq')
    ax.set_title('dose = '+ str(dose1) + 'um')
    ax.grid(True)
    #save figures
    plt.savefig(outf, bbox_inches=0)
    #if nMd < .25: #make graph colored with low rank score
	    #ax.patch.set_facecolor('#F1DDDD') #light red color
	    #plt.savefig(outf, facecolor=fig.get_facecolor(), edgecolor='white', bbox_inches=0)
    #else:
	    #plt.savefig(outf, bbox_inches=0)
    #fig.clear()

#loop through each unique perturbation - graph query ranks with dose
pearSet = []
spearSet = []
pertSet = set(qPert)
for pert in pertSet:
    iP = all_indices(pert, qPert)
    if len(iP) < 2:
	print pert + ' has only one instance'
	continue
    uDose = [qDose[i] for i in iP]
    uMetric = [medRnk[i] for i in iP] #metric drawn from large list
    fDose = [float(x) for x in uDose] #convert strings to float
    aDose = numpy.asarray(fDose) #convert to numpy array
    iD = aDose.argsort()
    sDose = [fDose[j] for j in iD] #sort local doses
    sMetric = [uMetric[j] for j in iD] #sort local metric
    Pear = scipy.stats.pearsonr([1,2,3,4],sMetric)
    Spear = scipy.stats.spearmanr([1,2,3,4],sMetric)
    pearSet.append(Pear)
    spearSet.append(Spear)
    ## plot 1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([1,2,3,4],sMetric,'o-')
    ax.set_ylim(0,1)
    ax.set_xticks([1,2,3,4])
    ax.set_xticklabels( tuple(sDose) )
    ax.set_ylabel('median normalized rank')
    ax.set_xlabel('dose')
    ax.set_title(pert + ' affogato self-connection rank')
    fname = pert + '_normalized_rank.png'
    outf = os.path.join(work_dir,fname)
    plt.savefig(outf, bbox_inches=0)
    
#return all indices where the input string matches the item in the list
def all_indices(value, qlist):
	indices = []
	indx = -1
	while True:
		try:
			indx = qlist.index(value, indx+1)
			indices.append(indx)
		except ValueError:
			break
	return indices


#scratch 

#graph alternative w/ pyplot
#f = pl.plot([1,2,3,4],uMetric,'ro')
#plt.setp(f, 'markersize', 10)
#plt.ylim(0,1)