import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd

work_dir = '/xchip/cogs/hogstrom/analysis/summly/bioAs'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

#retrieve sig ids from mongo - on per cp per cell line?
### make target_dict
# targetSheetF = '/xchip/cogs/projects/target_id/4June20    13/Informer2_drug_targets.txt'
# targetDict = {}
# pDescDict = {}
# with open(targetSheetF,'rt') as f:
#     for string in f:
#         splt = string.split('\r')
#         for i,line in enumerate(splt):
#             splt2 = line.split('\t')
#             pID = splt2[0] #the pert_id listed the line
#             pDesc = splt2[1]
#             targets = splt2[2]
#             targets = targets.split(';')
#             targets = [x for x in targets if x != '']
#             if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
#                 continue
#             else:
#                 targetDict[pID] = targets
#                 pDescDict[pID] = pDesc

cgsCells = ['A375', 'A549', 'ASC', 'HA1E', 'HEPG2', 'HCC515', 'HT29', 'MCF7', 'NPC', 'PC3', 'VCAP']
pertList = ['BRD-K68202742','BRD-K81418486','BRD-K12184916','BRD-K59369769']
for pert in pertList:
    outDir = work_dir + '/' + pert
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    file1 = outDir + '/' + 'sigs.grp'
    sigs = []
    for cell in cgsCells:
        CM = mu.CMapMongo()
        #dose min
        # pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell,'pert_dose':{'$gt':5}},{'sig_id':True},limit=1)
        #no dose min
        pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell},{'sig_id':True},limit=1)
        if pert_query:
            sigs.append(pert_query[0])
    with open(file1, 'w') as f:
        [f.write(x + '\n') for x in sigs]
    metric = 'wtcs'
    queryDir = outDir + '/sig_query'
    if not os.path.exists(queryDir):
        os.mkdir(queryDir)
    cmd = ' '.join(['rum -q local -f sig_query_tool',
             '--sig_id ' + file1,
             '--metric ' + metric,
             '--column_space full',
             '--out ' + queryDir,
             '--mkdir false',             
             '--save_tail false'])
             # '--row_space bing', 
    os.system(cmd)


trich_SigIds = '/xchip/cogs/hogstrom/analysis/summly/trichostatin-a/sig_ids.grp'
metric = 'wtcs'
cmd = ' '.join(['rum -q local -f sig_query_tool',
         '--sig_id ' + trich_SigIds,
         '--metric ' + metric,
         '--column_space full',
         '--out ' + work_dir,
         '--save_tail false'])
os.system(cmd)

#run summly as a system call
# matlab -r "A=3;MyScript"
# matlab -nodesktop -r "brd = 'BRD-K68065987'; /xchip/cogs/hogstrom/scripts/jailbird/summly/run_summly.m"
# matlab -nodesktop -r "/xchip/cogs/hogstrom/scripts/jailbird/summly/run_summly.m"
# matlab -r "run_summly(brd = 'BRD-K81418486')"
# to_run = '/xchip/cogs/hogstrom/scripts/jailbird/summly/run_summly.m'
# DosCmd = 'matlab -wait -automation -nosplash -r "run \'' + to_run + "', exit\""
# matlab -wait -automation -nosplash -r "brd = 'BRD-K68065987'; run /xchip/cogs/hogstrom/scripts/jailbird/summly/run_summly.m"

### run external queries
work_dir = '/xchip/cogs/hogstrom/analysis/summly/cflix_pannel_minus_phenox'
metric = 'wtcs'
fup = '/xchip/cogs/hogstrom/analysis/summly/cflix_pannel/cflix_up_minus_phenox.gmt'
fdn = '/xchip/cogs/hogstrom/analysis/summly/cflix_pannel/cflix_dn_minus_phenox.gmt'
cmd = ' '.join(['rum -q local -f sig_query_tool',
         '--uptag ' + fup,
         '--dntag ' + fdn,         
         '--metric ' + metric,
         '--column_space full',
         '--row_space full',
         '--out ' + work_dir,
         '--save_tail false'])
os.system(cmd)

### read in summly results

inFile = '/xchip/cogs/hogstrom/analysis/summly/bioAs/BRD-K12184916/summly_result/jul25/my_analysis.sig_summly_tool.2013072508560215/BRD-K12184916/BRD-K12184916_summly.txt'
sumRes = pd.io.parsers.read_csv(inFile,sep='\t')

# filter to only cps
pd.io.parsers.read_csv
cpRes = sumRes[sumRes['pert_type'] == 'trt_cp']
cpRes['rank'] = np.arange(1,len(cpRes)+1)
cgsRes = sumRes[sumRes['pert_type'] == 'trt_sh.cgs']
oeRes = sumRes[sumRes['pert_type'] == 'trt_sh.oe']

drugFile = '/xchip/cogs/projects/target_id/ctd2_annots/ctd2_merged_mapped_genes.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
targetSet = set(drugLabels['gene_dw_annot'])
bigGroups = ['AKT','BRAF','EGFR','HDAC','MAPK','MTOR','HSP90','PIK3','PPARG','TP53']

brdGrpList = []
for grp in bigGroups:
    grpPerts = drugLabels['pert_id'][drugLabels['gene_dw_annot'] == grp]
    brdGrpList.extend(grpPerts.values)
    for brd in grpPerts:
        indSum = cpRes[cpRes['pert_id'] == brd]['sum_score']
        sumScore = indSum.values[0]
        indrank = cpRes[cpRes['pert_id'] == brd]['rank']
        rank = indrank.values[0]
        percRnk = rank / float(len(cpRes))
        print str(percRnk) + ' ' + str(sumScore)

mtors = list(grpPerts.values)
cpRes[cpRes['pert_id'] == mtors[0]]
cpRes[cpRes['pert_id'] == mtors[0]]['pert_id']






