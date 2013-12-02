#! /usr/bin/env python
'''
load in the Rajiv's pcl sig_introspect results

clean up duplicate groups and duplicate BRDs

'''
import getpass
import pymongo
import numpy as np, pandas as pd
from cmap.util import debug
import cmap
from cmap.io import gmt
from os import path
import pandas as pd
import codecs
import cmap.io.gct as gct

from pymongo import MongoClient
from cmap.util.mongo_utils import CredentialsObject
import cmap.util.mongo_utils as mu

#introspect input
filePCLgrps = '/xchip/cogs/projects/pharm_class/pcl_shared_target.txt'
pclFrm = pd.io.parsers.read_csv(filePCLgrps,sep='\t')

#introspect results
introspectFile = '/xchip/cogs/web/icmap/summly_introspect/pcl_target_wtcs.combined/self_connectivity.txt'
iSpectFrm = pd.io.parsers.read_csv(introspectFile,sep='\t')
iSpectFrm = iSpectFrm.sort('median_rankpt',ascending=False)

mtrxFile = '/xchip/cogs/web/icmap/summly_introspect/pcl_target_wtcs.combined/self_rankpt_n758x758.gctx'
gt = gct.GCT(mtrxFile)
gt.read()
rnkptMtrx = gt.frame

#proportion of dup BRDs
iSpectFrm['unique_iname_proportion'] = np.nan
inameSet = set()
inameDict = {}
#loop through each groups
for icol in iSpectFrm.index:
    col = iSpectFrm.ix[icol]
    inameStr = col['pert_iname']
    inameList = inameStr.split('|')
    inameSet = inameSet.union(set(inameList)) #rolling set of pert_inames
    inameDict[col['group_id']] = set(inameList)
    LenIname = len(inameList)
    LenUnqIname = len(set(inameList))
    propUnique = (LenUnqIname)/float(LenIname) # proportion of inames that are unique
    iSpectFrm.ix[icol,'unique_iname_proportion'] = propUnique
outF = '/xchip/cogs/hogstrom/analysis/scratch/PCL_proportion_of_unique_inames.txt'
iSpectFrm.to_csv(outF,sep='\t',index=True,header=True)

#find compounds which overlap in multiple 
nInm = len(inameDict)
overlapCount = pd.DataFrame(np.zeros([nInm, nInm]),
    index=inameDict.keys(),
    columns=inameDict.keys())
overlapProp = pd.DataFrame(np.zeros([nInm, nInm]),
    index=inameDict.keys(),
    columns=inameDict.keys())
for iname1 in inameDict:
    grp1Set = inameDict[iname1]
    for iname2 in inameDict:
        grp2Set = inameDict[iname2]
        interSect = grp1Set.intersection(grp2Set)
        nIntersect = len(interSect)
        overlapCount.ix[iname1,iname2] = nIntersect #overlap counts
        propIntersect = nIntersect/float(len(grp1Set)) #proportion of overlap
        overlapProp.ix[iname1,iname2] = propIntersect
outF = '/xchip/cogs/hogstrom/analysis/scratch/PCL_overlap_proportion_matrix.txt'
overlapProp.to_csv(outF,sep='\t',index=True,header=True)

## what are the top overlaps between groups?
# grUpperMtrx = grtr_value_upper_mtrx(overlapProp.values)
# upperFrm = pd.DataFrame(grUpperMtrx,
#         index=overlapProp.index,
#         columns=overlapProp.index)
# overlapSer = upperFrm.unstack()
upperFrm = overlapProp.copy()
np.fill_diagonal(upperFrm.values, np.nan)
# overDiag = np.fill_diagonal(overlapProp.values, np.nan)
# upperFrm = pd.DataFrame(overDiag,
#         index=overlapProp.index,
#         columns=overlapProp.index)
overlapSer = upperFrm.unstack()
overlapSer = overlapSer[~overlapSer.isnull()] #remove nulls 
overlapSer.sort(ascending=False)
outF = '/xchip/cogs/hogstrom/analysis/scratch/PCL_overlap_proportion_list.txt'
overlapSer.to_csv(outF,sep='\t',index=True,header=True)
# overlapSer[overlapSer > .5]

#look at group overlap
# what proprtion of the group A overlap with group B
pclGrped = pclFrm.groupby('class')
pclCounts = pclGrped.size()
pclCounts[pclCounts > 3]


#first pass at reducing group size:
#1) if a PCL has <50% unique compounds - eliminate
iUniquePCLs = iSpectFrm['unique_iname_proportion'] > .5
uniquePCLs = iSpectFrm['group_id'][iUniquePCLs]
#2) if one PCL has >60% of its compounds in another PCL - choose the larger one
overlapThresh = .6
overThresh = overlapSer[overlapSer>overlapThresh]
removeList = []
for i in overThresh.index:
    iRev = (i[1],i[0])
    overlap1 = overlapSer[i]
    overlap2 = overlapSer[iRev]
    if overlap2 > overlap1:
        removeList.append(i[0])
removeSet = set(removeList)
# take out overlapping PCLs
uniqueSet = set(uniquePCLs.values)
pclKeep = uniqueSet - removeSet
# reindex, leaving groups with many brd dupes
# and compound redundancy
iSpectFrm.index = iSpectFrm['group_id']
g = list(pclKeep)
# g = ['GNAS','ADCY5']
reducedISpec = iSpectFrm.reindex(g)
medRnkpt = reducedISpec['median_rankpt'].copy()
medRnkpt.sort('median_rankpt', ascending=False)
# keepNames = medRnkpt.index.values
outF = '/xchip/cogs/hogstrom/analysis/scratch/pcl_keepers.txt'
medRnkpt.to_csv(outF,index=True,header=True,sep='\t')
curratedFile = '/xchip/cogs/hogstrom/analysis/scratch/pcl_keepers_currated.txt'
#use Rajiv's matlab validvar tool to change list strings
#write using mkgrp
nameModFile = '/xchip/cogs/hogstrom/analysis/scratch/pcl_keepers_mod_currated.txt'
nameMod = pd.io.parsers.read_csv(nameModFile,sep='\t',header=None)
nameMod[0].values

#write new index files with groups ordered by median rnkpt:
fileNames = []
for pcl in nameMod[0].values:
    fName = 'heatmap_' +pcl + '.png'
    fileNames.append(fName)    
indexfile = '/xchip/cogs/hogstrom/analysis/scratch/index_pcl_lh.html'
with open(indexfile,'w') as f:
    lineWrite = '<h2>PCLs by median rnkpt </h2>'
    f.write(lineWrite + '\n')
    for fName in fileNames:
        lineWrite =  '<img src=' + fName + '>'
        f.write(lineWrite + '\n')
        lineWrite2=  '<h6>' + fName + '</h6>'
        f.write(lineWrite2 + '\n')
def grtr_value_upper_mtrx(mtrx):
    '''
    -take an unsymetric matrix
    -put the larger of two values in an upper matrix
    
    '''        
    nm = len(mtrx)
    grtrMtrx = np.zeros((nm,nm))
    for i1 in range(nm):
        for i2 in range(nm):
            val1 = mtrx[i1,i2]
            val2 = mtrx[i2,i1]
            if val1 > val2:
                grtrMtrx[i1,i2] = val1
            else:
                grtrMtrx[i1,i2] = val2
    # grtrMtrxUp = np.triu(grtrMtrx,k=1)
    iUp = np.tril_indices(nm)
    grtrMtrx[iUp] = np.nan
    return grtrMtrx


onz = np.ones([4,4])
iUp = np.tril_indices(4, k=-2)
onz[iUp] = np.nan

np.fill_diagonal(onz, np.nan)