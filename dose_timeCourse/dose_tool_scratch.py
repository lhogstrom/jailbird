#! /usr/bin/env python

import cmap.analytics.sc as sc
import os.path
import cmap.analytics.signature_strength as ss
import numpy
import scipy
import cmap.io.gct as gct
#import ljh_dose_analysis_tool as dose

#plot tool
import pylab as pl
import matplotlib.pyplot as plt

cellLine = 'MCF7'
timeP = '24H'
gctfile = '/xchip/obelix/pod/brew/pc/ASG001_%s_%s/by_pert_id_pert_dose/ASG001_%s_%s_COMPZ.MODZ_SCORE_LM_n85x978.gctx' % (cellLine,timeP,cellLine,timeP) 

#make a gct object
db = gct.GCT()
db.read(gctfile)

#make sc object for signature strength
sco = sc.SC()
sco.add_sc_from_gctx_meta(gctfile)
ss = sco.s

#make signature for each dose
work_dir = '/xchip/cogs/hogstrom/analysis/scratch'
#work_dir = os.getcwd() #set work_dir var as pwd
fup = '/xchip/cogs/hogstrom/analysis/scratch/tmp_up_list.gmt'
fdn = '/xchip/cogs/hogstrom/analysis/scratch/tmp_dn_list.gmt'
open(fup,'w') #overwrite existing grp file
open(fdn, 'w') #overwrite existing grp file
n_edge = 50
db = gct.GCT()
db.read(gctfile)
cids = db.get_cids()
pertIDs = [x.split(':')[1] for x in cids]
doses = [float(x.split(':')[2]) for x in cids]
perts = db.get_column_meta('pert_desc')
probes = db.get_rids()
mtrx = db.matrix #matrix of data from gct file
#loop through each column of data
for i,pertID in enumerate(pertIDs):
    profile = mtrx[:,i]
    n_prof = len(profile)
    iprofile = profile.argsort() #indices that sort array
    iprofile = iprofile[::-1] #switch indicies to decend
    sprofile = profile[iprofile]
    itop = iprofile[0:(n_edge)]
    ibot = iprofile[-n_edge:n_prof]
    col_name = perts[i] + '_' + str(doses[i]) + 'um_' + cellLine + '_' + timeP
    ptop = [] 
    pbot = []
    for j,it in enumerate(itop):
	    ptop.append(probes[it]) #make probe id list
    for j,ip in enumerate(ibot):
	    pbot.append(probes[ip]) #make probe id list
    #to gmt list 
    #with open(os.path.join(work_dir,'up_list.grp'),'w') as f:
    #with open(os.path.join(work_dir,'tmp_up_lm_list.grp'),'w') as f:
    with open(fup,'a') as f:
	    f.write(col_name + '\t' + col_name + '\t')
	    for pt in ptop:
	        f.write(pt + '\t')
	    f.write('\n')
    with open(fdn,'w') as f:
	    f.write(col_name + '\t' + col_name + '\t')
	    for pb in pbot:
	        f.write(pb + '\t')
#python system call
cd work_dir
cmd = 'rum -q local query_tool --uptag ' + fup + ' --dntag ' + fdn + ' --metric eslm'
os.system(cmd)


## load in results
frslt = '/xchip/cogs/hogstrom/analysis/scratch/Nov20/dose_analysis_tool.1353449771597/nov20/my_analysis.query_tool.2012112017162991/result_ESLM.COMBINED_n85x398050.gctx'
rslt = gct.GCT()
rslt.read(frslt)

fname1 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_desc__affogato.txt';
affPert = open(fname1).read().splitlines()
fname2 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_id_affogato.txt';
affPertID = open(fname2).read().splitlines()
fname3 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/sig_ids_affogato.txt';
affSigID = open(fname3).read().splitlines()
rsltSigID = rslt.get_rids()

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

qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qDose = db.get_column_meta('pert_dose')
ESmat = rslt.matrix[irsltSigID,:] #order matrix acording to sigID sort
#iES = ESmat.argsort(axis=0) #sort ascending 
iES = ESmat.argsort(axis=0)[::-1] #sort ascending 
n_inst = len(iES[:,1])
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
    ax.set_title('normalized rank score = '+ str(nMd))
    ax.grid(True)
    if nMd < .25: #make graph colored with low rank score
	    ax.patch.set_facecolor('#F1DDDD') #light red color
	    plt.savefig(outf, facecolor=fig.get_facecolor(), edgecolor='white', bbox_inches=0)
    else:
	    plt.savefig(outf, bbox_inches=0)

#plot - colored graph
fname = cmpd1 + '_query_rank.png'
outf = os.path.join(work_dir,fname)
fig = plt.figure(figsize=(8.0, 2.0))
ax = fig.add_subplot(111)
# the histogram of the data
n, bins, patches = ax.hist(i2, 30, facecolor='green', alpha=0.75)
ax.set_xlabel('query rank')
ax.set_ylabel('freq')
ax.grid(True)
ax.patch.set_facecolor('#F1DDDD') #light red color
#plt.show()
plt.savefig(outf, facecolor=fig.get_facecolor(), edgecolor='white', bbox_inches=0)



#short plot code - median self connection ranks across compounds
fig = plt.figure()
ax = fig.add_subplot(111)
# the histogram of the data
n, bins, patches = ax.hist(medRnk, 20, normed=1, facecolor='green', alpha=0.75)
ax.set_xlabel('average normalized rank')
ax.set_ylabel('freq')
ax.grid(True)
#plt.show()
plt.savefig(outf, bbox_inches=0)

#plot ranks of self-connections
fig = plt.figure()
ax = fig.add_subplot(111)

# the histogram of the data
n, bins, patches = ax.hist(i2, 50, normed=1, facecolor='green', alpha=0.75)

# hist uses np.histogram under the hood to create 'n' and 'bins'.
# np.histogram returns the bin edges, so there will be 50 probability
# density values in n, 51 bin edges in bins and 50 patches.  To get
# everything lined up, we'll compute the bin centers
bincenters = 0.5*(bins[1:]+bins[:-1])
# add a 'best fit' line for the normal PDF
y = mlab.normpdf( bincenters, mu, sigma)
l = ax.plot(bincenters, y, 'r--', linewidth=1)

ax.set_xlabel('rank')
ax.set_ylabel('freq')
#ax.set_title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#ax.set_xlim(40, 160)
#ax.set_ylim(0, 0.03)
ax.grid(True)
plt.show()



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
	
	
	#sort sigIDs between query result and affogato saved
	iaffSigID = sorted(range(len(affSigID)), key = affSigID.__getitem__) #index of sorting
	irsltSigID = sorted(range(len(rsltSigID)), key = rsltSigID.__getitem__) #index of sorting
	
	
# run cmap pymongo query
import cmap.util.mongo_utils as mu
CM = mu.CMapMongo()
CM.find({'pert_desc':'AKT1'},{'sig_id':True},limit=10)
sigs = CM.find({'pert_desc':'AKT1'},{'sig_id':True},limit=10)
mu.sig_info(sigs,out='/xchip/cogs/cflynn/tmp/sig_info_example.txt')

CM = mu.CMapMongo()
qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qStr = qPertID[1]
sqStr = qStr[0:13]
cmpdSigIds = CM.find({'pert_id':sqStr},{'sig_id':True,'pert_id':True}) #list sig_id and pert_id of query result


