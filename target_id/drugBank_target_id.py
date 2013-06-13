import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo

work_dir = '/xchip/cogs/projects/target_id/DrugBank_12June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)
### make target_dict
# targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_short.txt'
# '/xchip/cogs/projects/target_id/4June2013/Informer2_drug_targets.txt'
targetSheetF = '/xchip/cogs/projects/target_id/7June2014/A2_DrugBank_targets_tab.txt'
targetDict = {}
pDescDict = {}
with open(targetSheetF,'rt') as f:
	for string in f:
		splt = string.split('\r')
		for i,line in enumerate(splt):
			splt2 = line.split('\t')
			pID = splt2[0] #the pert_id listed the line
			pDesc = splt2[1]
			targets = splt2[2:]
			targets = [x for x in targets if x != '']
			# targets = targets.split(';')
			if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
				continue
			else:
				targetDict[pID] = targets
				pDescDict[pID] = pDesc

### use pert_info collection to get sig_ids in mongo
# cellList = []
# for pert in targetDict.keys()[:10]:
# 	pertdb = mutil.CMapMongo(mongo_location = None, collection = 'pert_info')
# 	p1 = pertdb.find({'pert_id':'BRD-M79902621'},{'sig_id':True})
# 	g = p1[0]
# 	gSplit = g.split('\'')
# 	sigIDs = [x for x in gSplit if len(x) >= 5]
# 	cells = [x.split('_')[1] for x in sigIDs]
# 	cellList.extend(cells)
# print p1
# type(p1[0])

### which targets have CGS signatures
#get all CGS gene IDs
CM = mu.CMapMongo()
# pert_List = CM.find({'pert_type':{'$regex':pert}},{'sig_id':True,'cell_id':True})
CGSbyCell = CM.find({'pert_type':'trt_sh.cgs'},{'pert_iname':True})
CGSgeneSyms = set(CGSbyCell)
#check overlap with DB targets
nestedTargets = targetDict.values()
DBtargets = [item for sublist in nestedTargets for item in sublist]
setDBtargets = set(DBtargets)
DBcgsOverlap = setDBtargets.intersection(CGSgeneSyms)

targetDictCGS = {}
for pert in targetDict:
	for gene in targetDict[pert]:
		if gene in DBcgsOverlap:
			if targetDictCGS.has_key(pert):
				targetDictCGS[pert].append(gene)
			else:
				targetDictCGS[pert] = [gene]
#list of targets with CGS:
OverlapTargets = [item for sublist in targetDictCGS.values() for item in sublist]

test1 = 'OEB001_A375_96H:BRDN0000399163:-666' #set random sig_id to initialize dgo object
test2 = 'OEB001_A375_96H:BRDN0000400484:-666'

## test OE
# reload(dgo)
# dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_OE_connection')
# dg.add_dictionary(targetDict=targetDictCGS)
# dg.get_sig_ids(genomic_pert='OE')
# # dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
# dg.make_result_frames(gp_type='OE')
# dg.test_known_OE_connections(pDescDict=pDescDict,gp_type='OE')
# dg.FDR_correction(pDescDict=pDescDict)

### test KD
reload(dgo)
dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_KD_connection')
dg.add_dictionary(targetDict=targetDict)
# dg.get_sig_ids(genomic_pert='KD')
# dg.run_drug_gene_query(max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD')
# dg.test_known_connections(gp_type='KD',pDescDict=pDescDict,make_graphs=False)
dg.test_unknown_rank_product(gp_type='KD')
dg.FDR_correction(pDescDict=pDescDict)


#make table of connections pass FDR
for conn in dg.connectionsPassFDR:
	drug = conn[:13]
	pert = pDescDict[drug]
	gp = conn.split('-')[2]
	print drug + ':' + pert + ':' + gp

### exploratory connections - rank product
