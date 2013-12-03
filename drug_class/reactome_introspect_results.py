#! /usr/bin/env python
'''
load in the Rajiv's pcl sig_introspect results for reactome results

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

#introspect results
introspectFile = '/xchip/cogs/web/icmap/summly_introspect/canonical_pathway_wtcs.combined/self_connectivity.txt'
iSpectFrm = pd.io.parsers.read_csv(introspectFile,sep='\t')
iSpectFrm = iSpectFrm.sort('median_rankpt',ascending=False)


#proportion of dup BRDs
iSpectFrm['unique_iname_proportion'] = np.nan
inameSet = set()
inameDict = {}
inameSer = pd.Series()
#loop through each groups
for icol in iSpectFrm.index:
    col = iSpectFrm.ix[icol]
    inameStr = col['pert_iname']
    inameList = inameStr.split('|')
    inameSet = inameSet.union(set(inameList)) #rolling set of pert_inames
    inameDict[col['group_id']] = list(set(inameList))
    LenIname = len(inameList)
    LenUnqIname = len(set(inameList))
    propUnique = (LenUnqIname)/float(LenIname) # proportion of inames that are unique
    iSpectFrm.ix[icol,'unique_iname_proportion'] = propUnique
    tmpSer = pd.Series(col['group_id'],index=inameList)
    inameSer = pd.concat([inameSer,tmpSer])
outF = '/xchip/cogs/hogstrom/analysis/scratch/reactome_proportion_of_unique_inames.txt'
iSpectFrm.to_csv(outF,sep='\t',index=True,header=True)

# how many groups goes a given gene fall into?
inameSer.name = 'pathway'
inameFrm = pd.DataFrame(inameSer)
inameFrm['gene'] = inameFrm.index
geneGrped = inameFrm.groupby('gene')
geneCounts = geneGrped.size()
geneCounts = geneCounts.order(ascending=False)
outF = '/xchip/cogs/hogstrom/analysis/scratch/reactome_gene_freq_in_pathways.txt'
geneCounts.to_csv(outF,sep='\t',index=True,header=True)

groupSer = iSpectFrm['group_id']
outF = '/xchip/cogs/hogstrom/analysis/scratch/reactome_group_list.txt'
groupSer.to_csv(outF,sep='\t',index=False,header=False)


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
outF = '/xchip/cogs/hogstrom/analysis/scratch/reactome_overlap_proportion_matrix.txt'
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
outF = '/xchip/cogs/hogstrom/analysis/scratch/reactome_overlap_proportion_list.txt'
overlapSer.to_csv(outF,sep='\t',index=True,header=True)
# overlapSer[overlapSer > .5]


nameModFile = '/xchip/cogs/hogstrom/analysis/scratch/currated_reactome_group_list.txt'
nameMod = pd.io.parsers.read_csv(nameModFile,sep='\t',header=None)
nameMod[0].values

#write new index files with groups ordered by median rnkpt:
fileNames = []
for pcl in nameMod[0].values:
    fName = 'heatmap_' +pcl.lower() + '.png'
    fileNames.append(fName)    
indexfile = '/xchip/cogs/hogstrom/analysis/scratch/index_reactome_lh.html'
with open(indexfile,'w') as f:
    lineWrite = '<h2>PCLs by median rnkpt </h2>'
    f.write(lineWrite + '\n')
    for fName in fileNames:
        lineWrite =  '<img src=' + fName + '>'
        f.write(lineWrite + '\n')
        lineWrite2=  '<h6>' + fName + '</h6>'
        f.write(lineWrite2 + '\n')
