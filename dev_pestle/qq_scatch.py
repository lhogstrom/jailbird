import cmap.analytics.sc as sc
import os.path
import cmap.analytics.signature_strength as ss
import numpy
import scipy
import cmap.io.gct as gct
import cmap.util.mongo_utils as mutil
import matplotlib.pyplot as plt

#return all indices where the input string matches the item in the list
def _all_indices(value, qlist):
	indices = []
	indx = -1
	while True:
		try:
			indx = qlist.index(value, indx+1)
			indices.append(indx)
		except ValueError:
			break
	return indices

cellLine = 'MCF7'
timeP = '24H'
gctfile = '/xchip/obelix/pod/brew/pc/ASG001_%s_%s/by_pert_id_pert_dose/ASG001_%s_%s_COMPZ.MODZ_SCORE_LM_n85x978.gctx' % (cellLine,timeP,cellLine,timeP) 

#make a gct object
db = gct.GCT()
db.read(gctfile)

qPert = db.get_column_meta('pert_desc')
qPertID = db.get_column_meta('pert_id')
qDose = db.get_column_meta('pert_dose')
probeIDs = db.get_row_meta('id')

#set null distirbution of z-scores (currently normal)
ESmat = db.matrix
#calculate null distribution
mu, sigma = 0, 1
s = numpy.random.normal(mu, sigma, len(ESmat[:,1]))
s.sort()

# pertSet = set(qPert)
pertSet = set(qPertID)
# for pert in pertSet:
# 	iP = _all_indices(pert, qPertID) #index of doses on plate
# 	if len(iP) < 2:
# 		print pert + ' has only one instance'
# 		continue
# 	uDose = [qDose[i] for i in iP]
# 	fDose = [float(x) for x in uDose] #convert strings to float
# 	aDose = numpy.asarray(fDose) #convert to numpy array
# 	iD = aDose.argsort() #local ordering
# 	sDose = [fDose[j] for j in iD] #sort local doses
# 	iPo =  [iP[i] for i in iD] #ordered index
# 	#sMat = ESmat[:,iPo]
# 	#sMat.sort(axis=0)
# 	#mongo query for each unique pertID
# 	qStr = qPertID[iPo[0]] #set pertID
# 	if len(qStr) >= 13:
# 		qStr = qStr[0:13] #shorten qPertID
# 	CM = mutil.CMapMongo()
# 	#cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True})
# 	edge50Lst = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True,'up50_lm':True,'dn50_lm':True,'cell_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
# 	if isinstance(edge50Lst, int): #if the query result is a '0' int (eg empty), then skip
# 		continue
# 	nInstances = len(edge50Lst) #number of instances in db
# 	#count number of times a probe is in the top/bottom 50 genes of an instance
# 	upProbeCnts = [0] * len(probeIDs)
# 	dnProbeCnts = [0] * len(probeIDs)
# 	for j,inst in enumerate(edge50Lst):
# 		up50 = edge50Lst[j]['up50_lm']
# 		dn50 = edge50Lst[j]['dn50_lm']
# 		#loop through every gene in the top and bottom list - where does it live on the rank list?
# 		for prb in up50:
# 			if prb in probeIDs:
# 				iPrb = probeIDs.index(prb)
# 				upProbeCnts[iPrb] = upProbeCnts[iPrb] +1
# 		for prb in dn50:
# 			if prb in probeIDs:
# 				iPrb = probeIDs.index(prb)
# 				dnProbeCnts[iPrb] = dnProbeCnts[iPrb] +1
# 	#loop through each dose
# 	for d in iPo:
# 	#count probe enrichment and plot
# 			# cmpd1 = qPert[d]
# 			cmpd1 = qPertID[d]
# 			dose1 = qDose[d]
# 			zLst = db.matrix[:,d]
# 			iLst = zLst.argsort() #sort z-scores and save index
# 			sLst = zLst[iLst]
# 			sUpProbeCnts = [upProbeCnts[l] for l in iLst] #sort probe counts acording to z-score
# 			sDnProbeCnts = [dnProbeCnts[l] for l in iLst]
# 			#mkrs = numpy.sqrt(sprobeCnts) # non linear scaling of marker points
# 			sUpProbeCnts = [float(l) for l in sUpProbeCnts] #convert to float
# 			sDnProbeCnts = [float(l) for l in sDnProbeCnts] #convert to float
# 			# upPercMkrs = numpy.divide(sUpProbeCnts,max(sUpProbeCnts)) #divide by max count to make for relative frequency
# 			# dnPercMkrs = numpy.divide(sDnProbeCnts,max(sDnProbeCnts))
# 			upPercMkrs = numpy.divide(sUpProbeCnts,nInstances) #divide by total instances to make for relative frequency
# 			dnPercMkrs = numpy.divide(sDnProbeCnts,nInstances)
# 			upMkrs = numpy.multiply(upPercMkrs,100)
# 			dnMkrs = numpy.multiply(dnPercMkrs,100)
# 			#dummy plot to create legend dots, r1 b1
# 			fig2 = plt.figure()
# 			ax2 = fig2.add_subplot(111)
# 			r1 = ax2.plot(0,0,'r.',markersize=100,alpha=.25)
# 			b1 = ax2.plot(0,0,'b.',markersize=100,alpha=.25)
# 			plt.close()
# 			#plot real figure
# 			fig = plt.figure()
# 			ax = fig.add_subplot(111)
# 			ax.plot(s,s,'b')
# 			for j,sl in enumerate(sLst):
# 				ax.plot(s[j],sl,'r.',markersize=upMkrs[j],alpha=.25)
# 				ax.plot(s[j],sl,'b.',markersize=dnMkrs[j],alpha=.25)
# 			ax.set_ylabel('observed z-score')
# 			ax.set_xlabel('expected z-score')
# 			legStrUp = 'probe in 100% of ' + str(nInstances) + ' UP instances'
# 			legStrDn = 'probe in 100% of ' + str(nInstances) + ' DN instances'
# 			plt.legend([r1, b1], [legStrUp, legStrDn], numpoints=1, loc=4)
# 			#plt.legdend([b1], ['probe in 100% of ' + str(nInstances) + 'instances' ], numpoints=1)
# 			ax.set_title(pert + ' dose = ' + dose1)
# 			fname = pert + '_' + dose1 + 'um_connection_qq.png'
# 			outf = os.path.join(work_dir,fname)
# 			print 'saving graph for ' + pert + '_' + dose1
# 			plt.savefig(outf, bbox_inches=0)



# #one example
# # stuff needed:
# # marker size (1 and 2)
# # max marker size
# # vector of z-scores
# # null distribution
# #legend string 1
# #legend string 2
# #xlabel
# #ylabel
# #title
# #outputfile


# ##inputs
nullDist = numpy.random.normal(mu, sigma, len(ESmat[:,1]))
nullDist.sort()
zVec = zLst[iLst]
y_label = 'observed z-score'
x_label = 'expected z-score'
legStr1 = 'probe in 100% of ' + str(nInstances) + ' UP instances'
legStr2 = 'probe in 100% of ' + str(nInstances) + ' DN instances'
graph_title = pert + ' dose = ' + dose1
# outf = os.path.join(work_dir,fname)
upMkrs = numpy.multiply(upPercMkrs,100)
dnMkrs = numpy.multiply(dnPercMkrs,100)

# #graphing
fig, ax = plt.subplots(1)
# ax.plot(np.zeros_like(nullDist),linestyle='--',color='k')
r1 = ax.plot(0,0,'r.',markersize=100,alpha=.25)
b1 = ax.plot(0,0,'b.',markersize=100,alpha=.25)
plt.close()
#plot real figure
# fig = plt.figure()
# ax = fig.add_subplot(111)
fig, ax = plt.subplots(1)
ax.plot(nullDist,nullDist,'b')
for j,sl in enumerate(zVec):
	ax.plot(nullDist[j],sl,'r.',markersize=upMkrs[j],alpha=.25)
	ax.plot(nullDist[j],sl,'b.',markersize=dnMkrs[j],alpha=.25)
ax.set_ylabel(y_label)
ax.set_xlabel(x_label)
plt.legend([r1, b1], [legStr1, legStr2], numpoints=1, loc=4)
ax.set_title(graph_title)


#run with new qq tool
pert = 'BRD-K71879491-001-17-6'
iP = _all_indices(pert, qPertID) #index of doses on plate
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

d= iPo[0]
cmpd1 = qPertID[d]
dose1 = qDose[d]
zLst = db.matrix[:,d]
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
upMkrs = numpy.multiply(upPercMkrs,100)
dnMkrs = numpy.multiply(dnPercMkrs,100)

##define qq function inputs and run
nullDist = numpy.random.normal(mu, sigma, len(ESmat[:,1]))
nullDist.sort()
zVec = zLst[iLst]
y_label = 'observed z-score'
x_label = 'expected z-score'
legend_string1 = 'probe in 100% of ' + str(nInstances) + ' UP instances'
legend_string2 = 'probe in 100% of ' + str(nInstances) + ' DN instances'
graph_title = pert + ' dose = ' + dose1
# outf = os.path.join(work_dir,fname)
marker_size1 = numpy.multiply(upPercMkrs,100)
marker_size2 = numpy.multiply(dnPercMkrs,100)

import cmap.plot.qq as qq
qq.qq_plot(zVec,marker_size1,marker_size2,nullDist,graph_title,y_label, x_label, legend_string1,legend_string2)


