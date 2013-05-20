#!/bin/python
#run file like this
#execfile('/xchip/cogs/hogstrom/scripts/dose_timeCourse/ts_query.py')

#cd /xchip/cogs/hogstrom/scripts/dose_timeCourse

import pymongo

conn =  pymongo.Connection('vitalstatistix')
affogato = conn['affogato']
affogato.authenticate('cmap_user', 'l1000')
signature = affogato['signature']

#experiment with a given gene - hairpin
	testGene = 'NOTCH2'
	#qCurs = signature.find({'pert_type':'trt_sh','pert_desc':testGene},{'cell_id':True,'pert_type':True,'sig_id':True,'pert_desc':True,'is_gold':True})
	
	#query a bunch of stuff
	qCurs = signature.find({'pert_type':'trt_sh','pert_desc':testGene,'is_gold':1.0},{'cell_id':True,'pert_type':True,'sig_id':True,'pert_desc':True,'is_gold':True,'pert_dose':True,'pert_time':True,'cell_id':True})
	qLis = []
	for x in qCurs: 
	qLis.append(x)
	
	sigs = [str(x['sig_id']) for x in qLis]
	isGlds = [str(x['is_gold']) for x in qLis]
	dose = [str(x['pert_dose']) for x in qLis]
	pTime = [str(x['pert_time']) for x in qLis]
	cellID= [str(x['cell_id']) for x in qLis]

#Experiment with a given cell type - compound
cell1 = 'A375'
#filter for compounds, cellID, and is_gold
qCurs = signature.find({'pert_type':'trt_cp','cell_id':cell1,'is_gold':1.0},{'cell_id':True,'pert_type':True,'sig_id':True,'pert_desc':True,'is_gold':True,'pert_dose':True,'pert_time':True,'cell_id':True})
qLis = []
for x in qCurs: 
    qLis.append(x)

sigs = [str(x['sig_id']) for x in qLis]
isGlds = [str(x['is_gold']) for x in qLis]
dose = [str(x['pert_dose']) for x in qLis]
pTimes = [str(x['pert_time']) for x in qLis]
cellID= [str(x['cell_id']) for x in qLis]
cmpds = [str(x['pert_desc']) for x in qLis]

#save index of 24h time points
	pTs = ['24.0']
	i24h= []
	for i,pTime in enumerate(pTimes):
		if pTime == pTs[0] and cmpds[i] != '-666.0': #skip signatures without pert_desc
			i24h.append(i)
	
	h24 = []
	cmp24h = []
	for v in i24h:
		h24.append(pTimes[v])
		cmp24h.append(cmpds[v])

#save index of 6h time points
	pTs = ['6.0']
	i6h = []
	for i,pTime in enumerate(pTimes):
		if pTime == pTs[0] and cmpds[i] != '-666.0': #skip signatures without pert_desc
			i6h.append(i)
	
	h6 = []
	cmp6h = []
	for v in i6h:
		h6.append(pTimes[v])
		cmp6h.append(cmpds[v])

# match compounds in the 6h and 24h list
	#loop through each 6h exposure, find matching compound in 24h list
imatch24 = []
imatch6 = []
cmatch24 = []
cmatch6 = []
for i in i24h:
	for j in i6h:
		if cmpds[i] == cmpds[j]:
			imatch24.append(i)
			imatch6.append(j)
			cmatch24.append(cmpds[i])
			cmatch6.append(cmpds[j])

matchZip = zip(imatch24,imatch6,cmatch24,cmatch6)

#pTs = ['24.0','6.0']
#for pT in pTs:
	#print(pT)
	
#contained = [x for x in pTime if x in pTs[0]]
#[([int(item1 == item2) for item2 in list2], [n for n, item2 in enumerate(list2) if item1 == item2]) for item1 in list1]

#list(i[0] == i[1] for i in zip(list1, list2))

# write an outfile of sig IDs
with open('test_sig_ids.grp', 'w') as outfile:
for sig in sigs:
	outfile.write(str(sig) + '\n')
