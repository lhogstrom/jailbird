'''
Example script to run the summly_null anlaysis module

Larson Hogstrom, 12/2013
'''

import cmap.analytics.summly_null as SN

import numpy as np
import os
import cmap.util.mongo_utils as mu
import matplotlib.pyplot as plt
import pandas as pd
import cmap.io.gct as gct
from scipy import stats

reload(SN)

wkdir = '/xchip/cogs/projects/DOS/bioactivity_summary_Dec2013'
sn = SN.SummNull(out=wkdir)
sn.load_dmso_summ_results()
