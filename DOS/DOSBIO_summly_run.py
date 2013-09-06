#! /usr/bin/env python
'''
analyze the DOSBIO plates - combine with data in cmap database to perform query

use DOS signatures generate queries of the CGS data (cell line specific results)
'''

import os
import cmap.io.gct as gct
import cmap.analytics.sc as sc
import glob as glob
import cmap.util.mongo_utils as mutil
import cmap.util.progress as progress
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import cmap.analytics.dgo as dgo

wkdir = '/xchip/cogs/sig_tools/sig_summly/dosbio/list'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

### get DOS BIO cps from mongo
CM = mutil.CMapMongo()
# pert_List = CM.find({'sig_id':{'$regex':'DOSBIO'},'pert_iname':{'$regex':'BRD'}},{'sig_id':True,'pert_id':True,'pert_iname':True})
pert_List = CM.find({'sig_id':{'$regex':'DOSBIO'},'pert_iname':{'$regex':'BRD'}},{'pert_id':True})
dosbioSet = set(pert_List)
# check to make sure the brds are DOS compounds and don't represent known compounds
inameDict = {}
for brd in dosbioSet:
    inames = CM.find({'pert_id':brd},{'pert_iname':True})
    inameSet = set(inames)
    inameDict[brd] = inameSet

### get DOS BIO cps from mongo
CM = mutil.CMapMongo()
pert_List = CM.find({'sig_id':{'$regex':'DOSBIO'},'pert_iname':{'$regex':'BRD'}},{'pert_id':True})
dosbioSet = set(pert_List)
coreCells = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
inames = CM.find({'pert_id':{'$in':list(dosbioSet)},'sig_id':{'$regex':'DOS'},'cell_id':{'$in':coreCells},'is_gold':True},
            {'pert_id':True,'sig_id':True,'pert_iname':True,'cell_id':True},
            toDataFrame=True)
#are there inames that don't start with 
nonBRDinames = [x for x in inames['pert_iname'] if x[:3] != 'BRD']
nonDOSplateNames = [x for x in inames['sig_id'] if x[:3] != 'DOS']
cellSet = set(inames['cell_id'])


#are there perts with a suspisious amount of sig_ids? why would some DOS have 30 signatures?
sigCounts = inames.groupby('pert_id').count()
set(sigCounts['pert_id'])

#write one sig_id for each pert/cell line to a file:
sigCellGroup = inames.groupby(['pert_id','cell_id'])
sigTable = sigCellGroup.last() 
#limit to perts with is_gold in 4 or more cell lines
sigSummly = pd.DataFrame()
for brd in dosbioSet:
    if brd in inames['pert_id'].values:
        tmpTbl = sigTable.ix[brd]
        if len(tmpTbl) > 3:
            sigSummly = pd.concat([sigSummly,tmpTbl])
sigs = sigSummly['sig_id']
sigsF = wkdir + '/DOSBIO_summly_sigs.grp'
sigs.to_csv(sigsF,index=False)



for inx in sigTable.index:
    nT = len(sigTable.ix[inx])
    if nT > 3:

# print dos bios pert_ids to a file for rajiv
dosbioSer = pd.Series(list(dosbioSet))
sigsF = wkdir + '/DOSBIO_summly_perts.grp'
dosbioSer.to_csv(sigsF,index=False)

# load in Rajiv's file
sumSpaceFile = '/xchip/cogs/sig_tools/sig_query/results/a2_internal_wtcs.lm/query_info_n73597.txt'
summSpace = pd.read_csv(sumSpaceFile,sep='\t')

#make smaller summspace Frame
sumFrm = summSpace.copy()
sumFrm.index = summSpace['pert_id']
#sort based on pert_id
for brd in dosbioSet:
    if brd in sumFrm.index:
        tmpTbl = sumFrm.ix[brd] #table for the one pert_id
        tmpGrp = tmpTbl.groupby('cell_id')
        tmpLast = tmpGrp.last() # limit to only on sig_id per cell line
        tmpSer = tmpLast['sig_id'] # series of sig_ids
        brdFile = wkdir + '/' + brd + '.grp'
        tmpSer.to_csv(brdFile,index=False) # write sigs for each brd to a file


#submit summly jobs to lsf
for brd in dosbioSet:
    if brd in sumFrm.index:
        summlyMtrx = '/xchip/cogs/sig_tools/sig_query/results/a2_internal_wtcs.lm/'
        outDir = '/xchip/cogs/sig_tools/sig_summly/dosbio/summly_out'
        querySpace = wkdir + '/' + brd + '.grp'
        cmd = ' '.join(['rum -q hour -x sig_summly_tool',
                 summlyMtrx,
                 '--query_space ' + querySpace,
                 '--group_query true',
                 '--out ' + outDir])
        os.system(cmd)






