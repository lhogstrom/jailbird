#! /usr/bin/env python
'''
connect to Josh's mongoDB collection for drug-gene relationships

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

from pymongo import MongoClient
from cmap.util.mongo_utils import CredentialsObject
import cmap.util.mongo_utils as mu

CO = CredentialsObject()
CO.get_credentials()
server='23.23.153.188:27017,107.22.174.155:27017'
c = MongoClient(host=server) #,document_class='dex'
db = c['dex']
db.authenticate(CO.user,CO.pw)
# db.database['sig_info']
collection = db['targets']
collection.find_one()

g = collection.find()
# g = collection.find({'gene': 'HDAC1'})
dictList = list(g)
rFrame = pd.DataFrame(dictList)
rFrame.index = rFrame['_id']
# sourceSet = set(rFrame['source'])


### make a sheet for input for rajiv's tool
# class, pert_id, pert_iname, class_size

#DrugBank source
DBankFrm = rFrame[rFrame['source'] == 'DrugBank']
outFrm = DBankFrm.ix[:,['gene','pert_id']]
outFrm.columns = ['class','pert_id']
#eliminate duplicate columns
outFrm = outFrm.drop_duplicates(['class','pert_id'])

#TTD source
ttdFrm = rFrame[rFrame['source'] == 'TTD']
outFrm = ttdFrm.ix[:,['Name','pert_id']]
outFrm.columns = ['class','pert_id']

#set pert_inames
pIds = outFrm['pert_id'].values
CM = mu.CMapMongo()
inameRes = CM.find({'pert_id':{'$in':list(set(pIds))}}, #, 
        {'pert_id':True,'pert_iname':True},
        toDataFrame=True)
inameRes.index = inameRes['pert_id']
inameDict = inameRes['pert_iname'].to_dict()
inameSer = pd.Series(inameDict)
inameMtch = inameSer.reindex(outFrm['pert_id'])
outFrm['pert_iname'] = inameMtch.values

# set class_size
classGrped = outFrm.groupby('class')
classCounts = classGrped.size()
countMtch = classCounts.reindex(outFrm['class'])
outFrm['class_size'] = countMtch.values
outFrm = outFrm.sort('class_size',ascending=False)
# outF = '/xchip/cogs/projects/pharm_class/DrugBank_class_list_lh.txt'
outF = '/xchip/cogs/projects/pharm_class/TTD_class_list_lh.txt'
outFrm.to_csv(outF,sep='\t',index=False)



inameGrped = outFrm.groupby('pert_iname')
inameCount = inameGrped.size()
inameCount.sort(ascending=False) # sort inames acording to freq
# remove brd dups

#count how many groups an iname falls into
classCountDict = {}
for grp in inameGrped:
    grpName = grp[0]
    classSer = grp[1]['class']
    classCount = len(set(classSer))
    classCountDict[grpName] = classCount
classCountSer = pd.Series(classCountDict)
classCountSer.sort(ascending=False)
outF = '/xchip/cogs/hogstrom/analysis/scratch/compound_PCL_frequency.txt'
classCountSer.to_csv(outF,sep='\t',header=True)

#load in Introspect table

