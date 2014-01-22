'''
given input compounds, which ones are in summly space?

'''
import pandas as pd
import cmap.io.gct as gct

inFile = '/xchip/cogs/hogstrom/analysis/cmap_lung_pilot/cmap_compounds_short_list.txt'
inF = pd.read_csv(inFile,sep='\t')

# load summly matrix
summMtrx = '/xchip/cogs/projects/connectivity/summly/matrices/matched_mrp4_n7147x7147.gctx'
gt = gct.GCT()
gt.read(summMtrx)
summFrm = gt.frame

isSumm = inF['Structure ID'].isin(summFrm.index)
inF['is_summly'] = isSumm

outFile = '/xchip/cogs/hogstrom/analysis/cmap_lung_pilot/cmap_lung_is_summly.txt'
inF.to_csv(outFile,sep='\t')