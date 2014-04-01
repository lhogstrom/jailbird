'''
-run NMF jobs programatically

Larson Hogstrom, 4/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.io.gct as gct
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.util.progress as update
from matplotlib import cm
import glob
import shutil
import subprocess

# wkdir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/NMF_benchmark_development'
# if not os.path.exists(wkdir):
#     os.mkdir(wkdir)

# directory of NMF result prefix and matrix dimentions
# dimDict = {'LINCS_core_c9_LM':'n4716x978',
# 'LINCS_core_c9_bing':'n4713x10638',
dimDict = { 'PC3_c20_LM':'n585x978',
'PC3_c20_INF':'n585x10638',
'MCF7_c20_INF':'n652x10638',
'MCF7_c20_LM':'n652x978',
'MCF7_c9_INF':'n652x10638',
'MCF7_c9_LM':'n652x978',
'PC3_c9_INF':'n585x10638',
'PC3_c9_LM':'n585x978'}    

#move input files from one directory to another
# wkdir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation'
# wkdir2 = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation2'
# for prefix in dimDict:
#     print prefix
#     dim = dimDict[prefix]
#     mvDir = wkdir2 + '/' + prefix
#     if not os.path.exists(mvDir):
#         os.mkdir(mvDir)
#     gfile = '/' + prefix + '/clique_compound_classes_' + dim + '.gct'
#     aFile = '/' +  prefix + '/clique_compound_classes.v2.txt'
#     shutil.copy(wkdir+gfile, wkdir2+gfile)
#     shutil.copy(wkdir+aFile, wkdir2+aFile)

#specifications for subprocess
processes = set()
max_processes = 9 
### run jobs
wkdir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation2'
for prefix in dimDict:
    print prefix
    dim = dimDict[prefix]
    arg1 = wkdir + '/' + prefix # working directory
    arg2 = 'clique_compound_classes_' + dim
    cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/NMF_code_no_viz.v1.R', # 
         arg1,
         arg2])
    # os.system(cmd)
    processes.add(subprocess.Popen(cmd,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes if p.poll() is not None)
