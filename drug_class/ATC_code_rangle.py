'''
munge class label files
'''
import numpy as np
import os
import cmap
import pandas as pd
import matplotlib.pyplot as plt

drugFile = '/xchip/cogs/projects/cp_annot/ATC_codes/ATC.XLS'
# atcLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
# atcL = pd.ExcelFile(drugFile,kind='xls')

codeFile='/xchip/cogs/projects/cp_annot/ATC_codes/drugbank_cmap_ATCcode_20130930.csv'

indexFile='/xchip/cogs/projects/cp_annot/ATC_codes/2013_ATC_Index.csv'
