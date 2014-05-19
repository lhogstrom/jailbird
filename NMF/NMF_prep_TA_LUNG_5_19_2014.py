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

wkdir = '/xchip/cogs/projects/NMF/TA_lung_OE_May_2014/TA_OE_qnorm'

file_modz = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_COMPZ.MODZ_SCORE_n13968x22268.gctx'
file_qnorm = '/xchip/cogs/web/icmap/custom/TA/tnwork/datasets/for_jun10/TA_JUN10_QNORM_n37799x978.gctx'

gt = gct.GCT(src=file_qnorm)
gt.read()
qnorm = gt.frame

#split columns
colSer = pd.Series(qnorm.columns)
colSer.name = 'sig_id'
colSplit = colSer.str.split('_')

# annotate acording to sig_id fields
colFrame = pd.DataFrame(colSer)
colFrame['plate'] = colSplit.apply(lambda x: x[0])
colFrame['cell_line'] = colSplit.apply(lambda x: x[1])
colFrame['tp'] = colSplit.apply(lambda x: x[2])
colFrame['rep'] = colSplit.apply(lambda x: x[3])
colFrame['well'] = colSplit.apply(lambda x: x[4])

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
    cellFrm = qnorm.ix[:,sigs.values]
    nGt = gct.GCT()
    nGt.build_from_DataFrame(cellFrm)
    outF = cellDir + '/' + cell + '_TA_JUN10_COMPZ.MODZ_SCORE'
    nGt.write(outF,mode='gctx')

# convert gctx to gct
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

nComponents = 20
dimDict = {'A375':'n2245x978',
'A549':'n5608x978', # 
'AALE':'n3356x978',
'H1299':'n2597x978',
'HA1E':'n5371x978',
'PC3':'n2246x978',
'SALE':'n3215x978'}

#specifications for subprocess
processes = set()
max_processes = 9 
### run jobs
for cell in cell_grped.groups.keys():
    print cell
    dim = dimDict[cell]
    arg1 = wkdir + '/' + cell # working directory
    arg2 = cell + '_TA_JUN10_COMPZ.MODZ_SCORE_' + dim
    cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/NMF_code_no_viz.v2.R', # 
         arg1,
         arg2])
    # os.system(cmd)
    processes.add(subprocess.Popen(cmd,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes if p.poll() is not None)
