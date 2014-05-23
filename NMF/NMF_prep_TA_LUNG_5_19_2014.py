'''
organize data for TA lung cancer project

'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
import subprocess
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import cmap.io.gmt as gmt
import cmap.analytics.NMF_benchmarks as nmfb

# wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_May_2014/TA_OE_qnorm'
wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_May_2014/TA_OE_ZSPCINF'
# wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_May_2014/TA_OE_ZSPC_LM'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

################
### load data ##
################

file_modz = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_COMPZ.MODZ_SCORE_n13968x22268.gctx'
file_qnorm = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_QNORM_n37799x978.gctx'
file_zspcinf = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_ZSPCINF_n37799x22268.gctx'

gt = gct.GCT(src=file_zspcinf)
gt.read()
ds = gt.frame

processesed_type = 'ZSPCINF' # 'COMPZ.MODZ_SCORE', 'ZSPC_LM'
### reduce to LM genes 
# gLM = gct.GCT()
# gLM.read_gctx_row_meta(src=file_qnorm) # load file with LM genes
# lm_probes = gLM.get_rids()
# ds = ds.ix[ds.index.isin(lm_probes),:]

#split columns
colSer = pd.Series(ds.columns)
colSer.name = 'sig_id'
colSplit = colSer.str.split('_')

# annotate acording to sig_id fields
colFrame = pd.DataFrame(colSer)
colFrame['plate'] = colSplit.apply(lambda x: x[0])
colFrame['cell_line'] = colSplit.apply(lambda x: x[1])
colFrame['tp'] = colSplit.apply(lambda x: x[2])
colFrame['rep'] = colSplit.apply(lambda x: x[3])
colFrame['well'] = colSplit.apply(lambda x: x[4])

################################
### make cell line gct files ###
################################

# save matrix for each cell line in OE experiments
is_oe = colFrame.plate.str.match('TA.OE0')
oe = colFrame[is_oe]
cell_grped = oe.groupby('cell_line')
for grpT in cell_grped:
    cell = grpT[0]
    cellDir = wkdir + '/' + cell
    if not os.path.exists(cellDir):
        os.mkdir(cellDir)
    grp = grpT[1]
    sigs = grp.sig_id
    cellFrm = ds.ix[:,sigs.values]
    nGt = gct.GCT()
    nGt.build_from_DataFrame(cellFrm)
    outF = cellDir + '/' + cell + '_TA_JUN10_' + processesed_type
    nGt.write(outF,mode='gctx')

# convert gctx to gct
### run 'use Java-1.7' before running script ###
for cell in cell_grped.groups.keys():
    print(cell)
    cellDir = wkdir + '/' + cell
    outGCT = cellDir + '/' + cell
    globRes = glob.glob(outGCT+'*.gctx')
    print(globRes[0])
    cmd2 = 'convert-dataset -i ' + globRes[0]
    os.system(cmd2)

#########################
### Run NMF projection ##
#########################

# COMPZ.MODZ_SCORE
# nComponents = 20
# dimDict = {'A375':'n2245x978',
# 'A549':'n5608x978', # 
# 'AALE':'n3356x978',
# 'H1299':'n2597x978',
# 'HA1E':'n5371x978',
# 'PC3':'n2246x978',
# 'SALE':'n3215x978'}

# ZSPCINF
nComponents = 20
dimDict = {'A375':'n2245x22268',
'A549':'n5608x22268', # 
'AALE':'n3356x22268',
'H1299':'n2597x22268',
'HA1E':'n5371x22268',
'PC3':'n2246x22268',
'SALE':'n3215x22268'}

#specifications for subprocess
processes = set()
max_processes = 9 
### run jobs
for cell in cell_grped.groups.keys():
    print cell
    dim = dimDict[cell]
    arg1 = wkdir + '/' + cell # working directory
    # arg2 = cell + '_TA_JUN10_COMPZ.MODZ_SCORE_' + dim
    # arg2 = cell + '_TA_JUN10_ZSPCINF_' + dim
    arg2 = cell + '_TA_JUN10_' + processesed_type + '_' + dim
    cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/NMF_code_no_viz.v2.R', # 
         arg1,
         arg2])
    # os.system(cmd)
    processes.add(subprocess.Popen(cmd,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes if p.poll() is not None)

#########################
### load signature info ##
#########################

sFile = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/inst.info'
sigInfo = pd.read_csv(sFile,sep='\t')
sigInfo.index = sigInfo.distil_id

# important fileds in inst.info:
# pert_mfc_desc
# x_mutation_status
# x_preferredgenename
# x_tomconstructname

# reindex acording to OE plates
sigInfo = sigInfo.reindex(oe.sig_id)
# sigGrped = sigInfo.groupby(['cell_id','pert_mfc_desc'])
cellGrped = sigInfo.groupby('cell_id')
for cellTup in cellGrped:
    cell = cellTup[0]
    cellFrm = cellTup[1]
    cellDir = wkdir + '/' + cell
    outF = cellDir + '/OE_annotations.txt'
    # reformat sig_id
    cellFrm['mod_sig_id'] = cellFrm.distil_id.str.replace(':','.')
    cellFrm.index = cellFrm.mod_sig_id
    cellFrm.to_csv(outF,sep='\t')
    ### make gene signature groups - gmt file
    # geneGrped = cellFrm.groupby('pert_mfc_desc')
    geneGrped = cellFrm.groupby('x_mutation_status')
    gmtList = []
    for grp in geneGrped:
        gmtDictUp = {}
        gmtDictUp['id'] = grp[0]
        # gmtDictUp['desc'] = grp[0]
        gmtDictUp['desc'] = str(list(set(grp[1].x_mutation_status)))
        gmtDictUp['sig'] = list(grp[1].index.values)
        gmtList.append(gmtDictUp)
    gmtOut = cellDir + '/mutation_status_oe_sig_id.gmt'
    gmt.write(gmtList,gmtOut)

#########################
### run NMF benchmarks ##
#########################

for prefix in dimDict:
    dim = dimDict[prefix]
    path1 = wkdir + '/' + prefix
    prefix1 = prefix + '_TA_JUN10_'+ processesed_type + '_' + dim
    outdir = path1 + '/mutation_status_benchmark_graphs'
    source_dir = path1
    Hfile = prefix1 + '.H.k' + str(nComponents) + '.gct'
    # WFile = prefix1 + '.W.k' + str(nComponents) + '.gct'
    # MI_file_component = prefix1 + '.MI.k' + str(nComponents) + '.gct'
    # MI_file_inspace = prefix1 + '.MI.input_space.gct'
    anntFile = 'OE_annotations.txt'
    # groupFile = path1 + '/gene_oe_sig_id.gmt'
    groupFile = path1 + '/mutation_status_oe_sig_id.gmt'
    # run NMF module 
    reload(nmfb)
    self = nmfb.NMFresult(source_dir)
    self.set_output_dir(out=outdir)
    self.load_NMF_H_matrix(Hfile)
    self.load_annotations(anntFile,sig_col=0,drop_extra_signatures=True,signature_group_file=groupFile)
    self.group_component_maps(match_field='signatures')

