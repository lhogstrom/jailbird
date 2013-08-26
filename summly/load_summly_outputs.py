#! /usr/bin/env python
'''
load summly results into pyton
'''
import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd
import matplotlib


#regular summly table
inSum = '/xchip/cogs/hogstrom/analysis/summly/bioAs/BRD-K81418486/summly_result_saveFull/aug23/my_analysis.sig_summly_tool.2013082313052491/BRD-K81418486/BRD-K81418486_summly.txt'
sumRes = pd.io.parsers.read_csv(inSum,sep='\t')

#full summly table
inSumFull = '/xchip/cogs/hogstrom/analysis/summly/bioAs/BRD-K81418486/summly_result_saveFull/aug23/my_analysis.sig_summly_tool.2013082313052491/BRD-K81418486/BRD-K81418486_full.txt'
sumFull = pd.io.parsers.read_csv(inSumFull,sep='\t')

