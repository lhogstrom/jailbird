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

###############################
### retrieve screening data ###
###############################

# screening file 
sFile = '/xchip/cogs/projects/pharm_class/lhwork/HTS_results_for_L1000_compounds.csv'
hts = pd.read_csv(sFile)
# assayGrped = hts.groupby('ASSAY_NAME','ASSAY_OBS_ID')

####################
### assay counts ###
####################

# which screens are most common across drugs?
assayGrped = hts.groupby('ASSAY_NAME')
first = assayGrped.first()
# tested in what fraction of compounds?
assayCounts = assayGrped.apply(lambda x: len(set(x.pert_id)))
assayCounts.name = 'compounds_tested_in_assay'
aProportion = assayCounts/len(set(hts.pert_id)) # proportion of compounds tested in assay
aProportion.name = 'percent_compounds'
aFrm = pd.concat([assayCounts,aProportion],axis=1)
aFrm = aFrm.sort('percent_compounds',ascending=False)
outF = wkdir + '/HTS_assay_counts.txt'
aFrm.to_csv(outF)

##################################
### assay of apriori interest  ###
##################################

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

#################
### HTS marix ###
#################
# put HTS data in matrix format:
# 1) hit matrix: pert_id x assay_type
# 2) score matrix: pert_id x assay_type

assayGrped = hts.groupby(['ASSAY_NAME','ASSAY_OBS_ID'])
hitMatrix = pd.DataFrame()
scoreMatrix = pd.DataFrame()
for gname, group in assayGrped:
    gPgrp = group.groupby('pert_id')
    ## first value per compound
    # firstGroup = gPgrp.first() # take the first result for each compound
    # colName = firstGroup.ASSAY_DB_ID[0] + '-' + firstGroup.ASSAY_NAME[0] + '-' + str(firstGroup.ASSAY_OBS_ID[0])
    ## median values per compound
    firstGroup = gPgrp.median() # take the median result for each compound
    colName = group.ASSAY_DB_ID.values[0] + '-' + group.ASSAY_NAME.values[0] + '-' + str(group.ASSAY_OBS_ID.values[0])
    ### store hit status
    grpHits = firstGroup.hit
    grpHits = grpHits.replace(np.nan, 0) #replace nan with 0
    grpHits.name = colName
    grpHits = pd.DataFrame(grpHits)
    hitMatrix = pd.concat([hitMatrix,grpHits],axis=1)
    ### store score 
    grpScore = firstGroup.RESULT_VALUE
    grpScore = grpScore.replace(np.nan, 0) #replace nan with 0
    grpScore.name = colName
    grpScore = pd.DataFrame(grpScore)
    scoreMatrix = pd.concat([scoreMatrix,grpScore],axis=1)
outHit = wkdir + '/HTS_hit_matrix_median.txt'
hitMatrix.to_csv(outHit,sep='\t')
outScore = wkdir + '/HTS_score_matrix_median.txt'
scoreMatrix.to_csv(outScore,sep='\t')

################################
### cliqueselect specificity ###
################################

### sig_cliqueselect_tool specificity     
cFileFull = '/xchip/cogs/projects/DOS/PCL_comparison_Feb142014/feb14/my_analysis.sig_cliquescore_tool.2014021414592391/clique_median_n80x7147.gctx'
gt2 = gct.GCT()
gt2.read(cFileFull)
cliqFull = gt2.frame
cliqFull.index = gt2.get_row_meta('pert_id')
cliqNM = cliqFull.reindex(scoreMatrix.index)
#selectivity ratio
def specif_ratio(x):
    x = x.order(ascending=False)
    max1 = x[0]
    # max2 = x[1:5] 
    max2 = x[1] 
    topTwoRatio = max1/max2
    pHigher = (max1-max2)/max2
    #ratio:
    # highest value
    # median of next 5 highest?
    return pHigher
ratioSer = cliqNM.apply(specif_ratio,axis=1)
ratioSer.name = 's-score'
ratioSer.index.name = 'pert_id'
outF = wkdir + '/specificity_ratio.txt'
ratioSer.to_csv(outF,sep='\t',header=True)

######################
### HTS clustering ###
######################

### how similar are the compounds acording the the HTS 'sensors'?

mtrx = hitMatrix
Y = scipy.cluster.hierarchy.linkage(mtrx, method='centroid')
Z = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
cOrder = Z['leaves']
iPCL = mtrx.index[cOrder]
clustered = mtrx.reindex(index=iPCL)
#columns 
# Y_ = scipy.cluster.hierarchy.linkage(mtrx.T, method='centroid')
# Z_ = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
# cOrder_ = Z['leaves']
# iPCL_columns = mtrx.columns[cOrder]
# clustered = mtrx.reindex(index=iPCL,columns=iPCL_columns)

# make heatmap
plt.close()
# ccol.set_color_map()
fig = plt.figure(1, figsize=(20, 30))
plt.imshow(clustered.T.values,
    interpolation='nearest',
    aspect='auto')
xtickRange = range(0,clustered.shape[0])
xtcks = [x for x in clustered.index]
# ytickRange = range(0,clustered.shape[1])
# ytcks = [x for x in clustered.columns]
plt.xticks(xtickRange, xtcks,rotation=90)
# plt.yticks(ytickRange, ytcks)
plt.xlabel('compounds')
plt.ylabel('HTS assays')
# plt.yticks(np.arange(len(xtcks)),xtcks)
# plt.title('median connection of DOS to compound class')
plt.title('median connection of clique members to cliques')
plt.colorbar()
# outF = os.path.join(wkdir, 'dos_clique_heatmap_clust.png')
outF = os.path.join(wkdir, 'HTS_clustered_hit_heatmap.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

