'''
Find data related to doxyrubicin

5/15/2014
'''
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.io.gct as gct
import glob
import shutil
import cmap.util.mongo_utils as mu
import matplotlib.pyplot as plt
import numpy as np

### get doxyrubicin signatures
mc = mu.MongoContainer()
## KD
doxFrm = mc.sig_info.find({'pert_iname':'doxorubicin','pert_type':'trt_cp'},
            {'dn100_full':False,'up100_full':False,'dn100_bing':False,'up100_bing':False,'dn50_lm':False,'up50_lm':False,'pert_idose':False,'pert_dose_unit':False},toDataFrame=True)
## OE
cellGrped = doxFrm.groupby('cell_id')

# dose list for each cell line
cellGrped.apply(lambda x: set(x.pert_dose))
