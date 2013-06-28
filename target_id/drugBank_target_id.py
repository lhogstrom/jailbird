import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo

work_dir = '/xchip/cogs/projects/target_id/DrugBank_26June2013'
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


### test KD - two sided
reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_connection')
dg.add_dictionary(targetDict=targetDictCGS)
# dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')
# dg.test_known_connections(gp_type='KD',metric='spearman',pDescDict=pDescDict,make_graphs=True)
# dg.FDR_correction(pDescDict=pDescDict,metric='spearman',outName='apriori_connections_pass_FDR',alpha=0.2,make_graphs=False)
dg.test_known_connections(gp_type='KD',
                        metric='spearman',
                        pDescDict=pDescDict,
                        outName='two_sided_dg_graphs',
                        make_graphs=False,
                        n_rand=1000000,
                        connection_test='two_sided')
dg.FDR_correction(pDescDict=pDescDict,
                gp_type='KD',
                metric='spearman',
                outName='apriori_two_sided_pass_FDR',
                alpha=0.2,
                make_graphs=True)
dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR')
# # dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # # dg.test_unknown_rank_product(gp_type='KD')
# # # # dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')
outF = os.path.join(dg.outputdir,'drug-target_summary_two_sided.txt')
dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')
dg.store_parameters_rpt()

# ### test KD - one sided
# reload(dgo)
# dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_connection')
# dg.add_dictionary(targetDict=targetDict)
# # dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# # dg.run_drug_gene_query(metric='spearman',max_processes=10)
# # #wait until queries finish
# dg.make_result_frames(gp_type='KD',metric='spearman')
# # dg.test_known_connections(gp_type='KD',metric='spearman',pDescDict=pDescDict,make_graphs=True)
# # dg.FDR_correction(pDescDict=pDescDict,metric='spearman',outName='apriori_connections_pass_FDR',alpha=0.2,make_graphs=False)
# dg.test_known_connections(gp_type='KD',
#                         metric='spearman',
#                         pDescDict=pDescDict,
#                         outName='two_sided_dg_graphs',
#                         make_graphs=False,
#                         n_rand=1000000,
#                         connection_test='two_sided')
# dg.FDR_correction(pDescDict=pDescDict,
#                 gp_type='KD',
#                 metric='spearman',
#                 outName='apriori_two_sided_pass_FDR',
#                 alpha=0.2,
#                 make_graphs=True)
# dg.fdr_html_summary(fdrDir='apriori_two_sided_pass_FDR')
# # # dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # # # dg.test_unknown_rank_product(gp_type='KD')
# # # # # dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')
# outF = os.path.join(dg.outputdir,'drug-target_summary.txt')
# dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')
# dg.store_parameters_rpt()

# ### test KD - old
# # reload(dgo)
# dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_connection')
# dg.add_dictionary(targetDict=targetDict)
# # dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# # dg.run_drug_gene_query(metric='wtcs',max_processes=10)
# # #wait until queries finish
# dg.make_result_frames(gp_type='KD',metric='wtcs')
# dg.test_known_connections(gp_type='KD',metric='wtcs',pDescDict=pDescDict,make_graphs=False)
# dg.FDR_correction(pDescDict=pDescDict,metric='wtcs',outName='apriori_connections_pass_FDR',alpha=0.2)
# dg.fdr_html_summary(fdrDir='apriori_connections_pass_FDR')
# # dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='wtcs',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # dg.test_unknown_rank_product(gp_type='KD')
# # # dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')

# ### test KD - wtcs
# # reload(dgo)
# # dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_wtcs_is_gold')
# # dg.add_dictionary(targetDict=targetDict)
# # dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# # dg.run_drug_gene_query(metric='wtcs',max_processes=10)
# # #wait until queries finish
# dg.make_result_frames(gp_type='KD',metric='wtcs')
# dg.test_known_connections(gp_type='KD',metric='wtcs',pDescDict=pDescDict)
# dg.FDR_correction(pDescDict=pDescDict,metric='wtcs',outName='apriori_connections_pass_FDR',alpha=0.2)
# dg.fdr_html_summary(fdrDir='apriori_connections_pass_FDR')
# # dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='wtcs',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # dg.test_unknown_rank_product(gp_type='KD')
# # # dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')

# ### test KD - spearman
# # reload(dgo)
# # dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman_is_gold')
# # dg.add_dictionary(targetDict=targetDict)
# # dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# # dg.run_drug_gene_query(metric='spearman',max_processes=10)
# # #wait until queries finish
# dg.make_result_frames(gp_type='KD',metric='spearman')
# dg.test_known_connections(gp_type='KD',metric='spearman',pDescDict=pDescDict)
# dg.FDR_correction(pDescDict=pDescDict,metric='spearman',outName='apriori_connections_pass_FDR',alpha=0.2)
# dg.fdr_html_summary(fdrDir='apriori_connections_pass_FDR')
# # dg.gene_to_drug_similarity(testGene='ABCB5',gp_type='KD',metric='spearman',outName='gene_to_drug_connections',pDescDict=pDescDict,n_rand=10000,n_uncorrected=20)
# # # dg.test_unknown_rank_product(gp_type='KD')
# # # dg.FDR_correction(pDescDict=pDescDict,outName='FDR_pass_unknown')


# ### TEST OE
# # dg = dgo.QueryTargetAnalysis(test1,test2,work_dir + '/drug_OE_connection')
# # dg.add_dictionary(targetDict=targetDict)
# # dg.get_sig_ids(genomic_pert='OE')
# # # dg.run_drug_gene_query(max_processes=10)
# # # #wait until queries finish
# # dg.make_result_frames(gp_type='OE')
# # dg.test_known_OE_connections(pDescDict=pDescDict,gp_type='OE')
# # dg.FDR_correction(pDescDict=pDescDict,gp_type='OE',outName='FDR_pass_graphs_alpha2',alpha=0.2)

# # # #re-asign varabile so computation is not lost when reloading the class
# # # #dgCopy = dg
# reload(dgo)
# dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_wtcs_is_gold')
# dg.add_dictionary(targetDict=targetDict)
# dg.dfRank = dgCopy.dfRank
# dg.dfCS = dgCopy.dfCS
# dg.pVec = dgCopy.pVec
# dg.pDict = dgCopy.pDict
# dg.connectionsPassFDR = dgCopy.connectionsPassFDR

# #make table of connections pass FDR
# # for conn in dg.connectionsPassFDR:
# # 	drug = conn[:13]
# # 	pert = pDescDict[drug]
# # 	gp = conn.split('-')[2]
# # 	print drug + ':' + pert + ':' + gp

# ### exploratory connections - rank product
