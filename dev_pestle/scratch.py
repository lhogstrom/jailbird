
prog = progress.DeterminateProgressBar('Drug-target')
targetRanks = {}
targetRnkPercs = {}
targetCS = {}
for icell, cell1 in enumerate(cellList):
	celldir = os.path.join(work_dir,cell1) 
	outdir = os.path.join(work_dir,cell1,'sig_query_out')
	if not glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx'):
		print cell1 + 'no query result file'
		continue #if no results file, skip loop
	rsltF = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
	rslt = gct.GCT()
	rslt.read(rsltF)
	prog.update('analyzing {0}',icell,len(cellList))
	queryRids = rslt.get_rids()
	# for KD:
	# queryGenes = [x.split(':')[1] for x in queryRids]
	# for OE:
	queryGenes = [sigGeneDict[cell1][x] for x in queryRids]
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
	### find the indices of each query
	queryInd = {}
	for pert in pertSigIDs:
		iQuery = []
		iQuery = [i for i,x in enumerate(resPerts) if x[:13] == pert]
		if iQuery:
			queryInd[pert] = iQuery
	targetRanks[cell1] = {}
	targetRnkPercs[cell1] = {}
	targetCS[cell1] = {}
	sumF = os.path.join(celldir,cell1 + '_drug-target_connection_summary.txt')
	headers = ['query_sig','cp_pert_desc','target_KD_cgs','cs', 'query_rank', 'RnkPerc']
	with open(sumF,'w') as f:
		f.write('\t'.join(headers) + '\n')
		for pert in queryInd:
			targets = pertSigIDs[pert].keys()
			targetRanks[cell1][pert] = {}
			targetRnkPercs[cell1][pert] = {}
			targetCS[cell1][pert] = {}
			for target in targets:
				iTarget = [i for i,x in enumerate(queryGenes) if x == target]
				if iTarget:
					tarRanks = []
					tarCS = []
					tarRnkPercs = []
					for iTar in iTarget:
						iQuery = queryInd[pert]
						for iq in iQuery:
							cs = rslt.matrix[iTar,iq]
							rnk = rankMatrix[iTar,iq]
							RnkPerc = rnk/float(rslt.matrix.shape[0])
							tarCS.append(cs)
							tarRanks.append(rnk)
							tarRnkPercs.append(RnkPerc)
							#write summary table
							sig = resCids[iq]
							f.write('\t'.join([sig,pDescDict[pert],target,str(cs),str(rnk),str(RnkPerc)[:7]]) + '\n')
					targetCS[cell1][pert][target] = tarCS
					targetRnkPercs[cell1][pert][target] = tarRnkPercs
					targetRanks[cell1][pert][target] = tarRanks
	### put RnkPercs into a vector, test for FDR
	pVec = []
	for pert in targetRnkPercs[cell1]:
		RnkPercs = targetRnkPercs[cell1][pert].values()
		for RnkPerc in RnkPercs:
			pVec.extend(RnkPerc)
	#perform FDR on pVec
	[pBoolean, pCorrected] = fdr(pVec, alpha=0.1, method='indep')
	nPassFDR = sum(pBoolean)
	if nPassFDR:
		iPassFDR = [i for i,x in enumerate(pBoolean) if x == True]
		pPass = [pVec[i] for i in iPassFDR]
		pMaxThresh = max(pPass) #what is the largest RnkPerc that passed FDR
		#flag connections which pass FDR
		outF = os.path.join(celldir,cell1 + '_drug-target_FDR_summary.txt')
		headers = ['query_sig','cp_pert_desc','target_KD_cgs','cs', 'query_rank', 'RnkPerc']
		with open(outF,'w') as f:
			f.write('\t'.join(headers) + '\n')
			for pert in targetRnkPercs[cell1]:
				qInds = queryInd[pert]
				for target in targetRnkPercs[cell1][pert]:
					RnkPercs = targetRnkPercs[cell1][pert][target]
					for i,RnkPerc in enumerate(RnkPercs):
						if RnkPerc <= pMaxThresh:
							query = resCids[qInds[i]]
							cs = targetCS[cell1][pert][target][i]
							rank = targetRanks[cell1][pert][target][i]
							#query, target gene, cs, rank, RnkPerc
							f.write('\t'.join([query,pDescDict[pert],target,str(cs),str(rank),str(RnkPerc)[:7]]) + '\n')
	### check observed drug-target connections compared to null-distribution
	# of connections between drug and random gene knockdowns
	observedRanks = []
	observedCS = []
	for pert in targetRanks[cell1]:
		for target in targetRanks[cell1][pert]:
			observedRanks.extend(targetRanks[cell1][pert][target][:]) #create list of observed ranks
			observedCS.extend(targetCS[cell1][pert][target][:]) #create list of observed ranks
	n_rand = 1000
	randMtrxRank = np.zeros((len(observedRanks),n_rand))
	randMtrxCS = np.zeros((len(observedRanks),n_rand))
	for perm in range(n_rand):
		randRnkList = []
		randCSList = []
		for pert in queryInd:
			for target in targetRanks[cell1][pert]:
				nTargs = queryGenes.count(target)
				iRand = np.random.randint(rslt.matrix.shape[0], size=nTargs)
				for iTar in iRand:
					iQuery = queryInd[pert]
					cs = rslt.matrix[iTar,iQuery]
					rnk = rankMatrix[iTar,iQuery]
					randCSList.extend(cs)
					randRnkList.extend(rnk)
		randMtrxRank[:,perm] = randRnkList
		randMtrxCS[:,perm] = randCSList
	#plot ranks
	flat = randMtrxRank.flatten()
	h1 = plt.hist(flat,30,color='b',range=[0,max(flat)],label=['drug-random gene'],normed=True,alpha=0.5)
	h2 = plt.hist(observedRanks,30,color='r',range=[0,max(flat)],label='drug-target gene',normed=True,alpha=0.5)
	plt.legend()
	plt.xlabel('query rank')
	plt.ylabel('freq')
	plt.title('connection of drug and target gene knockdown')
	plt.savefig(os.path.join(celldir,cell1 + '_drug_target_rank_dist.png'))
	plt.close()
	#plot CS scores
	flat = randMtrxCS.flatten()
	h1 = plt.hist(flat,50,color='b',range=[min(flat),max(flat)],label=['drug-random gene'],normed=True,alpha=0.5)
	h2 = plt.hist(observedCS,50,color='r',range=[min(flat),max(flat)],label='drug-target gene',normed=True,alpha=0.5)
	plt.legend()
	plt.xlabel('cs score')
	plt.ylabel('freq')
	plt.title(cell1 + ' - connection of drug and target gene knockdown')
	plt.savefig(os.path.join(celldir,cell1 + '_drug_target_CS_dist.png'))
	plt.close()

### drug-gene pairs across all cell lines
pairDict = {}
for pert in pertSigIDs:
	if pert in noCMAPmatch['compound']:
		continue
	else:
		for target in targetDict[pert]:
			if target in noCMAPmatch['gene']:
				continue
			else:
				key1 = pert + ':' + target
				pairDict[key1] = []
				for cell1 in targetRnkPercs:
					if pert in targetRnkPercs[cell1].keys():
						if target in targetRnkPercs[cell1][pert]:
							percRanks = targetRnkPercs[cell1][pert][target]
							pairDict[key1].extend(percRanks)
				### graph 1% club - WaddenGram
				if pairDict[key1]: 
					if min(pairDict[key1]) < .01 or max(pairDict[key1]) > .99:
						print key1
						sKeysStr = []
						for i,cell1 in enumerate(targetRnkPercs):
							sKeysStr.append(cell1)
							if pert in targetRnkPercs[cell1].keys():
								if target in targetRnkPercs[cell1][pert]:
									percRanks = targetRnkPercs[cell1][pert][target]
									yVals = np.repeat(i+1,len(percRanks))
									plt.scatter(percRanks,yVals)
						plt.xlim((0, 1))
						plt.ylim((0,i+2))
						plt.yticks(range(1, i + 2), sKeysStr, rotation = 0)
						plt.xlabel('percent rank of drug in CGS query')
						plt.ylabel('cell line')
						plt.title(pDescDict[pert] + ' - ' + target)
						plt.savefig(os.path.join(work_dir,pert +'_' + target + '_connections.png'))
						plt.close()

### make html summary page
indexF = os.path.join(work_dir,'index.html')
with open(indexF,'w') as f:
	lineWrite = '<h2><CENTER>Connection between a drug and the knockdown of its gene target by cell line <CENTER></h2>'
	f.write(lineWrite + '\n')
	for cell1 in cellList:
		cellDir = os.path.join(work_dir,cell1)
		if os.path.exists(cellDir):
			lineWrite =  '<a href="' + cell1 + '/' + cell1 + '.html">' + cell1 + '</a> <BR>'
			f.write(lineWrite + '\n')
### make individual pages for each cell line
for cell1 in cellList:
	outdir = os.path.join(work_dir,cell1)
	if os.path.exists(outdir):
		pageF = os.path.join(outdir,cell1+ '.html')
		with open(pageF,'w') as f:
			csF = os.path.join(outdir,cell1 + '_drug_target_CS_dist.png')
			lineWrite = '<BR>' + '<tr><td><img src=' + csF + '></td></tr>'
			f.write(lineWrite + '\n')
			rankF = os.path.join(outdir,cell1 + '_drug_target_rank_dist.png')
			lineWrite = '<BR>' + '<tr><td><img src=' + rankF + '></td></tr>'
			f.write(lineWrite + '\n')
			#generate and write summary table
			lineWrite = '<h2><CENTER>Connections which pass FDR (based on query rank) <CENTER></h2>'
			f.write(lineWrite + '\n')
			sigF = os.path.join(outdir,cell1 + '_drug-target_FDR_summary.txt')
			if os.path.exists(sigF):
				table_data = []
				with open(sigF,'rt') as Fread:
					for line in Fread:
						table_data.append(line.split('\t'))
				htmlcodeFDR = HTML.table(table_data)
				f.write(htmlcodeFDR + '\n')
			else:
				f.write('<TD>No tested connections pass FDR</TD>')
			lineWrite = '<h2><CENTER>All connections tested <CENTER></h2>'
			f.write(lineWrite + '\n')
			### write connections tested
			sumF = os.path.join(outdir,cell1 + '_drug-target_connection_summary.txt')
			sumFrame = pd.read_table(sumF)
			sortSumFr = sumFrame.sort(columns='RnkPerc',ascending=True)
			sortSumFr = sortSumFr.set_index([range(len(sortSumFr))])
			table_data = []
			for i in range(len(sortSumFr)):
				line = [sortSumFr.query_sig[i],
						sortSumFr.cp_pert_desc[i],
						sortSumFr.target_KD_cgs[i],
						str(sortSumFr.cs[i]),
						str(sortSumFr.query_rank[i]),
						str(sortSumFr.RnkPerc[i])]
				table_data.append(line[:])
			htmlcode = HTML.table(table_data)
			f.write(htmlcode + '\n')
			### write connections not found in mongo
			lineWrite = '<h2><CENTER>Perturbations not found in CMAP database <CENTER></h2>'
			f.write(lineWrite + '\n')
			tblNoFind = []
			for pert in noCMAPmatch['compound']:
				tblNoFind.append([pert])
			lineWrite = '<h3>Compounds </h3>'
			f.write(lineWrite + '\n')
			htmlcode = HTML.table(tblNoFind)
			f.write(htmlcode + '\n')
			lineWrite = '<h3>collapsed gene signatures </h3>'
			f.write(lineWrite + '\n')
			tblNoFind = []
			for pert in noCMAPmatch['gene']:
				tblNoFind.append([pert])
			htmlcode = HTML.table(tblNoFind)
			f.write(htmlcode + '\n')


### scratch ### 

# # check rank results
# cids = rslt.get_cids()
# rids = rslt.get_rids()
# iquery = queryGenes.index('RARB')
# hpID = queryRids[iquery]
# iARF = rids.index(hpID)
# i743 = cids.index('CPC006_VCAP_24H:BRD-K93176058-001-02-8:10')
# rslt.matrix[iARF,i743]
# rankMatrix[iARF,i743]
# col1  = rslt.matrix[:,i743]
# ### get ranks to work correctly
# array = np.array([4,2,7,1])
# temp = array.argsort()[::-1]
# ranks = np.empty(len(array), int)
# ranks[temp] = np.arange(len(array)) + 1

# sortMatrix = np.argsort(rslt.matrix,axis=0)[::-1]
# rankMatrix = np.zeros_like(sortMatrix)
# for icol in range(rankMatrix.shape[1]):
# 	col = sortMatrix[:,icol]
# 	rankMatrix[col,icol] = np.arange(len(sortMatrix[:,icol])) + 1

# # rankMatrix = np.argsort(sortMatrix,axis=0) + 1
# rnkRslt = rnk.rnk_objects_from_gctx(str(rsltF))
# rankMatrix = rnkRslt.matrix
