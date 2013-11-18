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

c = MongoClient(host=server) #,document_class='dex'
db = c['dex']
db.authenticate(CO.user,CO.pw)
# db.database['sig_info']
collection = db['targets']
collection.find_one()


# mongo_location='current'
# class CredentialsObject(object):
#     '''class to store user credentials for CMap Mongo objects'''
#     def __init__(self,):
#         super(CredentialsObject, self).__init__()
    
#     def get_credentials(self):
#         '''
#         prompt the user for their password if access to CMap's internal 
#         authorization file is not allowed
#         '''
#         try:
#             mongo_servers = pd.read_table(path.join(cmap.mongo_path, 
#                                                     'mongo_servers.txt'),
#                                           index_col = 'mongo_location')
#             self.user = mongo_servers.loc['current', 'user_id']
#             self.pw = mongo_servers.loc['current', 'password']
#         except (IOError,ValueError):
#             self.user = raw_input('CMap Mongo User Name: ')
#             self.pw = getpass.getpass('CMap Mongo Password: ')

# CO = CredentialsObject()
# CO.get_credentials()

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


