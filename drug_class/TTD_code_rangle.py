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


# iFile = '/xchip/cogs/hogstrom/notes/TTD_annotations/from_CBNT_INCHIKEYs_added_unique_only.xlsx'
# iFile = '/xchip/cogs/hogstrom/notes/TTD_annotations/CBNT_INCHIKEYs.xlsx'
# cbnt = pd.ExcelFile(iFile,kind='xlsx')
# inchFrm = cbnt.parse('Sheet1')
iFile = '/xchip/cogs/hogstrom/notes/TTD_annotations/CBNT_INCHIKEYs.csv'
cbnt = pd.io.parsers.read_csv(iFile)
cbnt.index = cbnt['BROADID.1']

# ttdFrm = '/xchip/cogs/hogstrom/notes/TTD_annotations/ttd_cmap_export_20130925_inchikey3.xlsx'
# ttd = pd.ExcelFile(iFile,kind='xlsx')

ttdFrm = '/xchip/cogs/hogstrom/notes/TTD_annotations/TTD_targetID_pert_id_mapping_20130927.txt'
ttd = pd.io.parsers.read_csv(ttdFrm,sep='\t')

# specify directionality
brdLstLst = ttd.pert_ids_merged.str.split('|').tolist()
ttd['pert_ids_merged'] = brdLstLst
grouped = ttd.groupby(['TTD_target_ID','Category'])
ttd_cp_dict = {}
for group in grouped:
    groupName = group[1]['Name'].values[0]
    groupCat = group[1]['Category'].values[0]
    cpLstLst = group[1]['pert_ids_merged'].values
    cpLst = [item for sublist in cpLstLst for item in sublist] 
    cpLst = list(set(cpLst)) # make sure list is unique
    ttd_cp_dict[groupName+'-'+groupCat] = cpLst

## extend ttd to duplicate info for each pert_id
ttdE = pd.DataFrame()
lineList = []
for ix in ttd.index:
    lineSer = ttd.ix[ix]
    nBrds = lineSer['pert_ids_merged']
    if len(nBrds) == 1:
        brd = lineSer['pert_ids_merged']
        lineSer['pert_ids_merged'] = brd[0]
        lineList.append(lineSer)
    else:
        for brd in lineSer['pert_ids_merged']:
            lineSerOneBrd = lineSer.copy()
            lineSerOneBrd['pert_ids_merged'] = brd
            lineList.append(lineSerOneBrd)
ttdE = pd.DataFrame(lineList)
ttdE.index = ttdE['pert_ids_merged']

#make a dictionary for inames with more than one brd
inameGrped = ttdE.groupby('Drug_name')
GrpedDict = inameGrped.groups
inameDict = {}
multiDict = {} # which inames have more than one brd
for iname in GrpedDict:
    brdList = GrpedDict[iname]
    brdSet = list(set(brdList))
    inameDict[iname] = brdSet
    if len(brdSet) > 1:
        multiDict[iname] = brdSet

# for inames with multiple brds - chech inchy keys and names
dupBrdFrm = pd.DataFrame()
for ttd_iname in multiDict:
    brds = multiDict[ttd_iname]
    cbntRes = cbnt.ix[brds]
    cbntRes['ttd_iname'] = ttd_iname
    dupBrdFrm = pd.concat([dupBrdFrm,cbntRes],axis=0)
oFile = '/xchip/cogs/hogstrom/notes/TTD_annotations/ttd_brd_duplicates.csv'
dupBrdFrm.to_csv(oFile,sep='\t')


for iname in inameDict:
    ix = inameDict[iname]
    tmpFrm = ttdE.ix[ix]
    brdSet = list(set(tmpFrm['pert_ids_merged']))
    if len(brdSet):
        cbntNames = cbnt.ix[brdSet]    

