#! /usr/bin/env python

'''
-Analyze EMT signatures from Wai Leong Tam et. al. 2013
-run CMAP query on EMT signatures

LH 07/2014
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import cmap.util.mongo_utils as mu
import cmap.io.gct as gct

# load differential expression sigature from GEO2R
rFile = '/xchip/cogs/hogstrom/bathe/emt_signature/TWIST_vs_ctrl_vec.txt'
twst = pd.read_csv(rFile,sep='\t')


for x in range(0,200):
    print twst.ix[x,'Gene.symbol']


# run up/dn signatures through MSIGDB - do they capture other emt signatures? 
