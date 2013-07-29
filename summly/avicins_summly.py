import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd

work_dir = '/xchip/cogs/projects/avicins/gold_10um'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

avicinsBrds = ['BRD-A15100685','BRD-A33746814','BRD-A69592287','BRD-A70150975'] #avicin-d, avicin-g, oxetane, hydroxyl,
pDescDict = {'BRD-A15100685':'avicin-d','BRD-A33746814':'avicin-g','BRD-A69592287':'oxetane','BRD-A70150975':'hydroxyl'}

########################################
### count how many is gold there are for each compound - any dose
########################################

# for pert in avicinsBrds:
#     #dose min
#     sigAll = CM.find({'pert_id':{'$regex':pert},'pert_dose':{'$gt':3}},{'sig_id':True})
#     sigGold = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'pert_dose':{'$gt':3}},{'sig_id':True})
#     #no dose min
#     # sigAll = CM.find({'pert_id':{'$regex':pert}},{'sig_id':True})
#     # sigGold = CM.find({'pert_id':{'$regex':pert},'is_gold':True},{'sig_id':True})
#     nAll = len(sigAll)
#     nGold = len(sigGold)
#     percGold = nGold/float(nAll)
#     print ' '.join([pert,str(nAll),str(nGold), str(percGold)])

########################################
### write sig ids to file - run query ##
########################################
CM = mu.CMapMongo()
cellAll = CM.find({'pert_id':{'$regex':avicinsBrds[0]},'pert_dose':{'$gt':3}},{'cell_id':True})
cgsCells = list(set(cellAll))
# cgsCells = ['A375', 'A549', 'ASC', 'HA1E', 'HEPG2', 'HCC515', 'HT29', 'MCF7', 'NPC', 'PC3', 'VCAP', 'U266', 'JURKAT']
file1 = work_dir + '/' + 'avicin_sigs_10um.grp'
for pert in avicinsBrds:
    sigs = []
    for cell in cgsCells:
        #dose min - gold
        pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell,'pert_dose':{'$gt':9}},{'sig_id':True},limit=1)
        #dose min 
        # pert_query = CM.find({'pert_id':{'$regex':pert},'cell_id':cell,'pert_dose':{'$gt':9}},{'sig_id':True},limit=1)
        # no dose min
        # pert_query = CM.find({'pert_id':{'$regex':pert},'is_gold':True,'cell_id':cell},{'sig_id':True},limit=1)
        if pert_query:
            sigs.append(pert_query[0])
        # else:
        #     print 'no sig for ' + pert + ' in ' + cell
    with open(file1, 'a') as f:
        [f.write(x + '\n') for x in sigs]
### run query
metric = 'wtcs'
queryDir = work_dir + '/sig_query_10um'
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

# PI-3K/AKT,Stat-3,WNT/B-catenin, mTOR, NF-kB, cMYC, snail
# These pathways obviously support growth, proliferation, and/or survival
# glycolytic lines like Jurkat, PC-3, U266, and HA1E may have the most interesting profiles
# Metabolic enzymes like GAPDH or Aldolase or Actin isoforms like Arp2/3 
# surface in the query. Proteins with reactive cysteines are of great interest.


### read in summly results
graphDir = work_dir + '/graph_out'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
nCp = len(avicinsBrds)
avicinMtrx_sum_score = np.zeros((nCp,nCp))
avicinMtrx_PercSummly = np.zeros((nCp,nCp))
avicinMtrx_rank = np.zeros((nCp,nCp))
for ibrd,brd in enumerate(avicinsBrds):
    # brd = 'BRD-A15100685'
    basePath = work_dir + '/sig_query_10um'
    dateID = 'jul29/my_analysis.sig_summly_tool.2013072914520591'
    inFile = '/'.join([basePath,
                    dateID,
                    brd,
                    brd+'_summly.txt'])
    sumRes = pd.io.parsers.read_csv(inFile,sep='\t')
    # filter to only cps / cgs
    pd.io.parsers.read_csv
    cpRes = sumRes[sumRes['pert_type'] == 'trt_cp']
    cpRes['rank'] = np.arange(1,len(cpRes)+1)
    cgsRes = sumRes[sumRes['pert_type'] == 'trt_sh.cgs']
    cgsRes['rank'] = np.arange(1,len(cgsRes)+1)
    oeRes = sumRes[sumRes['pert_type'] == 'trt_sh.oe']
    ### matrix of avicin self connections
    for ibrd2, brd2 in enumerate(avicinsBrds):
        indSum = cpRes[cpRes['pert_id'] == brd2]['sum_score']
        if not indSum:
            print brd + ' ' + brd2 + ' not compared' 
            continue
        sumScore = indSum.values[0]
        indrank = cpRes[cpRes['pert_id'] == brd2]['rank']
        rank = indrank.values[0]
        percSummly = rank / float(len(cpRes))
        avicinMtrx_sum_score[ibrd,ibrd2] = sumScore
        avicinMtrx_PercSummly[ibrd,ibrd2] = percSummly
        avicinMtrx_rank[ibrd,ibrd2] = rank
        print brd + ' ' + brd2 + ' ' + str(percSummly) + ' ' + str(sumScore)
    # define genes of interest
    goi = ['GAPDH', 'ARP2', 'ARP3', 'ALDOA', 'ALDOB'] #genes of interes
    goi = ['PIK3CA', 'PIK3CB', 'AKT1', 'AKT2', 'MTOR', 'NFKB1','MYC', 'AMPK'] #genes of interes
    goiSet = set(goi)
    goiShort = goiSet.copy()
    # how many of these genes have been knocked down?
    cgsGenes = cgsRes['pert_iname'] 
    for gene in goi:
        if not gene in cgsGenes.values:
            print gene + ' not in CGS list'
            goiShort.remove(gene) #leave only genes we have profiled
    # place to find it
    for gene in goiShort:
        indSum = cgsRes[cgsRes['pert_iname'] == gene]['sum_score']
        sumScore = indSum.values[0]
        indrank = cgsRes[cgsRes['pert_iname'] == gene]['rank']
        rank = indrank.values[0]
        percSummly = rank / float(len(cgsRes))
        print gene + ' ' + str(percSummly) + ' ' + str(sumScore)
    #plot distribution of sum_scores
    metric = 'sum_score'
    metVals = sumRes[metric]
    metricVals = metVals.copy()
    metricVals.sort()
    metricVals = metricVals[::-1]
    plt.plot(metricVals)
    plt.xlabel('perturbagens')
    plt.ylabel(metric)
    plt.title('distributions of scores for ' + pDescDict[brd])
    plt.savefig(os.path.join(graphDir,brd + '_' + metric + '.png'))
    plt.close()
    # plt.show()

#make heatmap of rank score
plt.imshow(avicinMtrx_sum_score,interpolation='nearest')
plt.colorbar()
plt.show()
plt.close()

