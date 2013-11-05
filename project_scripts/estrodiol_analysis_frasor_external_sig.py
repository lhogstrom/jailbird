'''
analyze estrodiol compounds 

October 2013
'''
import numpy as np
from cmap.analytics.pert_explorer import PertExplorer
from cmap.analytics.cluster import HClust
import cmap.analytics.sc as sc
import cmap
import os
from os import path
from matplotlib import cm
import cmap.util.mongo_utils as mu
import subprocess
import cmap.tools.sig_dose_tool as sdt
import cmap.io.gct as gct
import pandas as pd
import cmap.io.gmt as gmt
from scipy import stats
import matplotlib.pyplot as plt

# get directory
dir1 = '/xchip/cogs/projects/pharm_class' 
wkdir = dir1 + '/estrodiol_analysis_Nov5'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

pDescDict = {'BRD-A60070924':'17-alpha-estradiol',
    'BRD-K04218075':'clomiphene',
    'BRD-K95309561':'dienestrol',
    'BRD-K45330754':'diethylstilbestrol',
    'BRD-K04046242':'equilin',
    'BRD-A74907996':'equol',
    'BRD-K18910433':'estradiol-17-beta',
    'BRD-K17016787':'estriol',
    'BRD-K81839095':'estrone',
    'BRD-A83237092':'fulvestrant',
    'BRD-K43797669':'genistein',
    'BRD-K63828191':'raloxifene',
    'BRD-K93754473':'tamoxifen',
    'BRD-K67174588':'toremifene'}
estrodiolBrds = pDescDict.keys()

### run sig_query_tool
frasUp = '/xchip/cogs/projects/pharm_class/estrodiol_analysis_Nov5/cflix_up_Frasor.gmt'
frasDn = '/xchip/cogs/projects/pharm_class/estrodiol_analysis_Nov5/cflix_dn_Frasor.gmt'
# run sig_query tool using the gene setls
# sig_query_tool('uptag', UP, 'dntag', DN, 'metric', METRIC)

# cmd = ' '.join(['rum -q local -f sig_query_tool',
#          '--sig_id ' + sigF,
#          '--metric ' + metric,
#          '--column_space custom',
#          '--cid ' + cidF,
#          '--out ' + outdir,
#          '--mkdir false',
#          '--save_tail false'])

outdir = wkdir + '/Frasor_sig_query_output'
if not os.path.exists(outdir):
    os.mkdir(outdir)
cmd = ' '.join(['rum -q hour -f sig_query_tool',
         '--uptag ' + frasUp,
         '--dntag ' + frasDn,
         '--out ' + outdir,
         '--row_space bing', #lm, bing, full
         '--metric wtcs',
         '--mkdir false',
         '--save_tail false'])

### run summly
summlyMtrx = wkdir + '/Frasor_sig_query_output'
outdir = wkdir + '/Frasor_sig_summly_output_bing'
if not os.path.exists(outdir):
cmd = ' '.join(['rum -q hour -x sig_summly_tool',
         summlyMtrx,
         '--query_space ' + querySpace,
         '--group_query true',
         '--out ' + outDir])
os.system(cmd)


