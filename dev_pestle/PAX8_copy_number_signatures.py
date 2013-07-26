import cmap.io.gct as gct
import glob, HTML
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection as fdr
import cmap.util.mongo_utils as mu
from cmap.tools import sig_slice_tool
from cmap.io import gct,plategrp,rnk
import cmap.util.progress as progress
import subprocess
import cmap.util.tool_ops as to
import cmap.analytics.dgo as dgo
# import cmap.analytics.oracle as oracle

work_dir = '/xchip/cogs/projects/target_id/PAX8_copy_number_15Jul'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

#write CTD2 sig ID files to disc
targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_drug_targets.txt'
targetDict = {}
pDescDict = {}
with open(targetSheetF,'rt') as f:
    for string in f:
        splt = string.split('\r')
        for i,line in enumerate(splt):
            splt2 = line.split('\t')
            pID = splt2[0] #the pert_id listed the line
            pDesc = splt2[1]
            targets = splt2[2]
            targets = targets.split(';')
            targets = [x for x in targets if x != '']
            if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
                continue
            else:
                targetDict[pID] = targets
                pDescDict[pID] = pDesc

### test KD
reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
dg.add_dictionary(targetDict=targetDict)
dg.get_sig_ids(genomic_pert='KD',is_gold=True)

# cnFile = '/xchip/cogs/projects/repurposing/PAX8/PAX8_50_50_s2n_n1x22268.gct'
# cnFile = '/xchip/cogs/projects/repurposing/PAX8/PAX8_100_100_s2n_n1x22268.gct'
cnFile = '/xchip/cogs/projects/repurposing/cn/cn_loss_lm_n1003x978.gctx'
cn = gct.GCT()
cn.read(cnFile)
#generate query - using sig_score file
metric='spearman'
max_processes=10
processes = set()
work_dir = dg.outputdir
for cell1 in dg.cell_lines_tested:
    cellDir = os.path.join(work_dir,cell1) 
    cidF = glob.glob(cellDir + '/' + cell1 + '_genomic_sig_ids_n*.grp')
    if not cidF:
        continue
    cidF = cidF[0]
    outdir = os.path.join(work_dir,cell1,'sig_query_out')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # sigF = os.path.join(cellDir, cell1 + '_genomic_sig_ids_n' + str(nCGS) + '.grp')
    sigF = os.path.join(cellDir,cell1 + '_cp_sig_ids.grp')
    # sig_query_tool('sig_score','score_s2n_n1x22277.gct', 'metric', 'wtcs', 'gset_size', 50,'row_space', 'full')
    cmd = ' '.join(['rum -q local -f sig_query_tool',
             '--sig_score ' + cnFile,
             '--metric ' + metric,
             '--column_space custom',
             '--cid ' + sigF,
             '--out ' + outdir,
             '--mkdir false',
             '--save_tail false'])
    processes.add(subprocess.Popen(cmd,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes if p.poll() is not None)

### make result frame 
# dg.make_result_frames(gp_type='OE',metric='spearman')
gp_type = 'KD'
work_dir = dg.outputdir
metric = 'spearman'
#which cell lines have a result dir
cellDirs = [f for f in os.listdir(work_dir) if os.path.isdir(work_dir+'/'+f)]
prog = progress.DeterminateProgressBar('dataframe read')
df = pd.DataFrame()
dfRank = pd.DataFrame()
#loop through each cell line add to df
for icell, cell1 in enumerate(cellDirs):
    #define directories and load in outputs
    outdir = os.path.join(work_dir,cell1,'sig_query_out')
    if not glob.glob(outdir + '/result_*.gctx'):
        print cell1 + ' no query result file'
        continue #if no results file, skip loop
    if metric == 'wtcs':
        rsltFile = glob.glob(outdir + '/result_WTCS.LM.COMBINED_n*.gctx')[0]
    if metric == 'spearman':
        rsltFile = glob.glob(outdir + '/result_SPEARMAN_n*.gctx')[0]
    rslt = gct.GCT()
    rslt.read(rsltFile)
    prog.update(cell1,icell,len(cellDirs))
    rsltF = rslt.frame
    # rsltF = rsltF.T
    indVals = rsltF.index.values
    pertVals = [ind.split(':')[1][:13] for ind in indVals]
    #make the column name gene and pert time
    geneVals = rsltF.columns
    # geneVals = []
    # for ind in rsltF.columns:
    #     if gp_type == 'KD':
    #         gene = ind.split(':')[1]
    #     if gp_type == 'OE':
    #         brdn = ind.split(':')[1]
    #         gene = dg.BRDNdict[brdn]
    #     tp = ind.split(':')[0].split('_')[-1]
    #     gname = '_'.join([gene, tp])
    #     geneVals.append(gname)
    if len(geneVals) > len(set(geneVals)):
        print 'duplicate CGS for this celline'
    newF = rsltF
    newF.index = [pertVals, rsltF.index.values]
    if gp_type == 'KD':
        newF.columns = geneVals
    if gp_type == 'OE':
        newF.columns = [geneVals, rsltF.columns.values]
    rankF = newF.rank(ascending=False,axis=1)
    perRankF = rankF / float(rankF.shape[1]) * 100.0
    #add cell line result to combined df
    if len(df) == 0:
        df = newF
        dfRank = perRankF
    else:
        df = pd.concat([df,newF],axis=0)
        dfRank = pd.concat([dfRank,perRankF],axis=0)
dg.dfCS = df
dg.dfRank = dfRank

dg.gene_to_drug_similarity(testGene='score',
                            gp_type='KD',
                            metric='spearman',
                            outName='/gene_to_drug',
                            pDescDict=pDescDict,
                            n_rand=10000,
                            n_uncorrected=20,
                            connection_test='two_sided')


# ### look at oncoDome genes
# inFile = '/xchip/cogs/projects/oncoDome/OncoDome_genes.txt'
# outDir = 'oncoDome_15July'
# if not os.path.exists(dg.outputdir+'/'+outDir):
#     os.mkdir(dg.outputdir+'/'+outDir)
# cgsList = []
# with open(inFile,'rt') as f:
#     for string in f:
#         splt = string[:-1]
#         cgsList.append(splt)
# geneAll = set(cgsList)
# # check to see gene has a CGS
# cnGenes = []
# for gene in geneAll:
#     if gene in dg.dfRank.columns:
#         cnGenes.append(gene)
# for gene in cnGenes:
#     dg.gene_to_drug_similarity(testGene=gene,
#                                 gp_type='KD',
#                                 metric='spearman',
#                                 outName=outDir + '/gene_to_drug',
#                                 pDescDict=pDescDict,
#                                 n_rand=10000,
#                                 n_uncorrected=20,
#                                 connection_test='two_sided')

# ### test known connections
# dg.test_known_connections(gp_type='KD',
#                         metric='spearman',
#                         pDescDict=pDescDict,
#                         outName='test_dg_graphs2',
#                         conn_thresh=.05,
#                         make_graphs=False,
#                         n_rand=100000,
#                         connection_test='two_sided')
# dg.FDR_correction(pDescDict=pDescDict,
#                 gp_type='KD',
#                 metric='spearman',
#                 outName='FDR_pass',
#                 alpha=0.2,
#                 make_graphs=True,
#                 specificity_graph=True)
# dg.fdr_html_summary(fdrDir='FDR_pass',specificity_graph=True)
# dg.store_parameters_rpt()
# outF = os.path.join(dg.outputdir,'drug-target_summary_peyton.txt')
# dg.make_target_summary(outF,dir_loc='apriori_connections_pass_FDR')
