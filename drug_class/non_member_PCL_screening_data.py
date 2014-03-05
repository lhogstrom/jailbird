'''
- examine high throughput screening data for cliq hits

Larson Hogstrom, 3/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.util.progress as update
from matplotlib import cm
import cmap.plot.colors as ccol
import scipy.cluster

wkdir = '/xchip/cogs/projects/pharm_class/lhwork/non_member_screening_data_3March2014'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

########################
### retrieve screening data ###
########################

# screening file 
sFile = '/xchip/cogs/projects/pharm_class/lhwork/HTS_results_for_L1000_compounds.csv'
hts = pd.read_csv(sFile)
# assayGrped = hts.groupby('ASSAY_NAME','ASSAY_OBS_ID')
assayGrped = hts.groupby('ASSAY_NAME')
first = assayGrped.first()

# 'HDAC2BroadScreen'
# 'A1 Apoptosis'

# which screens are most common across drugs?
# aCounts = assayGrped.count().PCL_hit_group
# aCounts = aCounts.order(ascending=False)
# tested in what fraction of compounds?
assayCounts = assayGrped.apply(lambda x: len(set(x.pert_id)))
assayCounts.name = 'compounds_tested_in_assay'
aProportion = assayCounts/len(set(hts.pert_id)) # proportion of compounds tested in assay
aProportion.name = 'percent_compounds'
aFrm = pd.concat([assayCounts,aProportion],axis=1)
aFrm = aFrm.sort('percent_compounds',ascending=False)
outF = wkdir + '/HTS_assay_counts.txt'
aFrm.to_csv(outF)

### apriori hits
# # 'HSP90_inhibitor',
hsp90_assays = ['HeatShockModulation.FluorProt.HSE-GFP', 'HeatShockModulation.FluorDye.Resazurin', 'MLPCNHeatShockFactor-1']
hspFrm = hts[hts.ASSAY_NAME.isin(hsp90_assays)]
hspHit = hspFrm[hspFrm.hit == 1]
# PI3K hits
pi3k_assays = ['PI3K/mTORModulation.FluorDye.Resazurin', 'PI3KinaseSignaling.Absorb.600']
pikFrm = hts[hts.ASSAY_NAME.isin(pi3k_assays)]
pikHit = pikFrm[pikFrm.hit == 1]
# 'HDAC_inhibitor',
hdac_assays = ['ChromatinBiology.ImmFluor.Eu-TRF', 'ChromatinBiology.FluorDye.FluorDeLys', 'HDAC2BroadScreen', 'MLPCNHDAC3']
hdacFrm = hts[hts.ASSAY_NAME.isin(hdac_assays)]
hdacHit = hdacFrm[hdacFrm.hit == 1]
# all apriori
apriori_assays = hsp90_assays + pi3k_assays + hdac_assays
aprioriFrm = hts[hts.ASSAY_NAME.isin(apriori_assays)]
aprioriHit = aprioriFrm[aprioriFrm.hit == 1]
aprioriHit = aprioriHit.sort('ASSAY_NAME')
outF = wkdir + '/apriori_HTS_assay_hits.txt'
aprioriHit.to_csv(outF,sep='\t',index=False)

### make a drug x assay type matrix
pertGrped = hts.groupby('pert_id')
for gname, group in pertGrped:
    g=1







