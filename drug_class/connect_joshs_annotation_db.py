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

sourceSet = set(rFrame['source'])
DBankFrm = rFrame[rFrame['source'] == 'DrugBank']

### make a sheet for input for rajiv's tool
# class, pert_id, pert_iname, class_size

#DrugBank source
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





#how many unique drug-gene pairs?
drugGenePairs = rFrame['pert_id'] + ':' + rFrame['gene']
setDGP = set(drugGenePairs)

geneGrped = DBankFrm.groupby('gene')

CatSer = rFrame['Category']
catFrm = rFrame[~pd.isnull(CatSer)]
np.isnan(rFrame['Category'].values)
rFrame['Category'] != np.nan



# ### copy fields from mongo_utils
#     _dflt_fields = {entry['id'] : entry['sig'] 
#                     for entry in gmt.read(path.join(cmap.mongo_path, 'default_fields.gmt'))}
#     # get mongo server information
#     _mongo_servers = pd.read_table(path.join(cmap.mongo_path, 
#                                              'mongo_servers.txt'),
#                                    index_col = 'mongo_location')

#     def __init__(self,CredentialsObject=CO, mongo_location='current', server='23.23.153.188:27017,107.22.174.155:27017',
#         database='a2y13q1',collection='sig_info'):
#         '''
#         constructor
#         '''
#         super(CMapMongo, self).__init__()
#         if mongo_location:
#             with open('/xchip/cogs/data/vdb/mongo/mongo_servers.txt','r') as f:
#                 headers = f.readline().strip().split('\t')
#                 mongo_location_idx = headers.index('mongo_location')
#                 database_idx = headers.index('sig_db')
#                 collection_idx = headers.index('sig_collection')
#                 replica_set_idx = headers.index('replica_set')
#                 lines = f.readlines()
#                 d = {}
#                 for line in lines:
#                     fields = line.strip().split('\t')
#                     mongo_loc = fields[mongo_location_idx]
#                     server = fields[replica_set_idx]
#                     database = fields[database_idx]
#                     collection = fields[collection_idx]
#                     d.update({mongo_loc: [server,database,collection]})

#             self.conn = pymongo.Connection(d[mongo_location][0])
#             self.database = self.conn[d[mongo_location][1]]
#             self.database.authenticate(CredentialsObject.user,CredentialsObject.pw)
#             self.collection = self.database[d[mongo_location][2]]

# server = '23.23.153.188:27017,107.22.174.155:27017'
# jdbConn = pymongo.Connection(d[mongo_location][0])
# jdbDatabase = jdbConn[d[mongo_location][1]]
# # jdbDatabase = jdbConn['dex']
# jdbDatabase.authenticate(CO.user,CO.pw)
# jdbCollection = jdbDatabase.database['sig_info']
# # jdbCollection = jdbDatabase.database['targets']\\

# # pymongo.mongo_client.MongoClient


