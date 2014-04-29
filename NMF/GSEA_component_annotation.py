'''
-Load in W matrix from NMF
-prep gene sets for GSEA

Larson Hogstrom, 4/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.util.progress as update

source_dir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c20_INF'
Wfile = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation/PC3_c20_INF/clique_compound_classes_n585x10638.W.k20.gct'
Wmtrx = pd.read_csv(Wfile,sep='\t',skiprows=[0,1],index_col=0) #,header=True
Wmtrx = Wmtrx.drop('Description',1)
# Wmtrx = Wmtrx.T

###
for nC1 in Wmtrx.columns:
    # nC1 = 'c11'
    c1 = Wmtrx.ix[:,nC1]
    c1 = c1.order()
    probSer = pd.Series(c1.index[-100:])
    outF = source_dir + '/PC3_' + nC1 + '_INF_n100_high_gene_weights.grp'
    probSer.to_csv(outF,index=False)
