'''
-Load in W matrix from NMF
-prep gene sets for GSEA

Larson Hogstrom, 4/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.util.progress as update
from matplotlib import cm

source_dir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c20_INF'
Wfile = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c20_INF/clique_compound_classes_n585x10638.W.k20.gct'
Wmtrx = pd.read_csv(Wfile,sep='\t',skiprows=[0,1],index_col=0) #,header=True
Wmtrx = Wmtrx.drop('Description',1)
# Wmtrx = Wmtrx.T

###
for nC1 in Wmtrx.columns:
    # nC1 = 'c11'
    c1 = Wmtrx.ix[:,nC1]
    c1 = c1.order()import cmap.io.gct as gct

    probSer = pd.Series(c1.index[-100:])
    outF = source_dir + '/PC3_' + nC1 + '_INF_n100_high_gene_weights.grp'
    probSer.to_csv(outF,index=False)


### plot of matrices
wkdir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c20_LM'
Wfile = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c20_LM/clique_compound_classes_n585x978.H.k20.gct'
# Wfile = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c20_LM/clique_compound_classes_n585x978.gct'
Wmtrx = pd.read_csv(Wfile,sep='\t',skiprows=[0,1],index_col=0) #,header=True
Wmtrx = Wmtrx.drop('Description',1)
# 
# fig = plt.figure(figsize=(20, 10), dpi=50)
plt.imshow(Wmtrx,
    interpolation='nearest',
    cmap=cm.gray_r,
    vmax=15)
outF = os.path.join(wkdir,'H_mtrx.png')
plt.savefig(outF, bbox_inches='tight')
plt.close()




