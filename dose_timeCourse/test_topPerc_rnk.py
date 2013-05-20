pertSet = set(qPert)
# pertSet = set(['DMSO', 'LY-294002', 'alpha-estradiol', 'alvespimycin'])
for pert in pertSet:
	cmpd1 = pert
	iP = _all_indices(pert, qPert) #index of doses on plate
	if len(iP) < 2:
		print pert + ' has only one instance on the plate'
		continue
	uDose = [qDose[i] for i in iP]
	fDose = [float(x) for x in uDose] #convert strings to float
	aDose = numpy.asarray(fDose) #convert to numpy array
	iD = aDose.argsort() #local ordering
	sDose = [fDose[j] for j in iD] #sort local doses
	iPo =  [iP[i] for i in iD] #ordered index
	qStr = qPertID[iPo[0]] #set pertID
	if len(qStr) >= 13:
		qStr = qStr[0:13] #shorten qPertID
	#run pymongo query
	CM = mutil.CMapMongo()
	#cmpdSigIds = CM.find({'pert_id':qStr},{'sig_id':True})
	cmpdSigIds = CM.find({'pert_id':{'$regex':qStr}},{'sig_id':True}) #search for the BRD-xxxxxxxxxxx within the pert_id field in the db
	n_cmpdInst = len(cmpdSigIds)
	if len(cmpdSigIds) < 1:
		print cmpd1 + ' has no instances in the cmap database'
		continue
	if len(cmpdSigIds) == 1:
		print cmpd1 + ' has one instance in the cmap database'
	#loop through each dose ()
	topTable = []
	for d in iPo:
	#count probe enrichment and plot
			cmpd1 = qPert[d]
			dose1 = qDose[d]
			iE = iES[:,d] #ES sort index for one column
			sSigID = []
			for y in iE:
				sSigID.append(rsltSigID[y]) #make sorted sig ID list
			i1 = [sSigID.index(y) for y in cmpdSigIds] #where instances of the compound of interest sit on the rank list
			i2 = numpy.array(i1) #convert list to numpy array
			avr = sum(i2)/len(i2) #what is the average ES rank
			md = numpy.median(i2) # what is the median ES rank
			nAv = float(avr)/n_inst #normalize acording to number of instances in db
			nMd = float(md)/len(iES[:,1]) #normalized median
			i1.sort()
			np = 1000 #define the top of the rank
			ntop = [x for x in i1 if x <= np]
			# PercTop = ntop/n_cmpdInst #percent of instances in the top of the rank
			nPr = float(len(ntop))/(len(i1)) #percent of instances at the top of the list
			topTable.append(nPr)
			prRnk.append(nPr)
			avRnk.append(nAv) #store average ES rank
			medRnk.append(nMd)
			#plotourse/test_topPerc_rnk.py
			# fname = cmpd1 + '_' + dose1 + '_query_rank.png'
			# outf = os.path.join(work_dir,fname)
			# fig = plt.figure(figsize=(8.0, 2.0))
			# ax = fig.add_subplot(111)
			# # the histogram of the data
			# n, bins, patches = ax.hist(i2, 30, facecolor='green', alpha=0.75)
			# #ax.set_xlim(0, n_inst)
			# ax.set_xlim(0, int(round(n_inst,-5))) #round instances to nearest 100k
			# ax.set_xlabel('query rank')
			# ax.set_ylabel('freq')
			# ax.set_title('dose = '+ str(dose1) +'um')
			# ax.grid(True)
			# plt.savefig(outf, bbox_inches=0)
	#top list for each compound 
	with open(qSumF,'a') as f:
		f.write(cmpd1 + '\t')
		for pr in topTable:
			f.write(str(pr) + '\t')
		f.write('\n')



## template 
gcto = db
cids = gcto.get_gctx_cid()
pert_descs = gcto.get_column_meta('pert_desc')
perts = [x.split(':')[1] for x in cids]
pert_desc_dict = dict(zip(perts,pert_descs))
unique_perts = list(set(perts))
unique_perts.sort()
unique_pDescs = list(pert_desc_dict.values())
unique_pDescs.sort()

# buld an environment for jinja2
cmap_base_dir = '/'.join(os.path.dirname(cmap.__file__).split('/')[0:-1])
env = jinja2.Environment(loader=jinja2.FileSystemLoader(cmap_base_dir + '/templates'))


index_page_template = env.get_template('Link_List_Template2_ljh.html')
# index_links = [pert_desc_dict[x] + '_detail.html' for x in unique_perts]
index_links = [x + '_detail.html' for x in unique_pDescs]
with open(os.path.join(work_dir,'index2.html'),'w') as f:
	f.write(index_page_template.render(title='Dose Analysis Results',
										links=zip(unique_pDescs,index_links),
										labels=unique_perts))
