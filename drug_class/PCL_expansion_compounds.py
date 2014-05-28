'''
-look up compound pubchem ids (CID) based on common name

Larson Hogstrom, 5/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import os

## load Pubchem synonyms and CIDs
gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/CID-Synonym-filtered'
cid = pd.read_csv(gFile,sep='\t',header=None,names=['CID','synonym'])

# load synonyms of interest:
sFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs.csv'
# sFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs.txt'
newPCL = pd.read_csv(sFile)
synonyms = newPCL.ix[:,"compound name"].values

################
### Get CIDS ### 
################

# find the synonymims that have CID matches
is_match = cid.synonym.isin(synonyms)
cid_match = cid[is_match]
mtch_grped = cid_match.groupby('synonym')
mtch_first = mtch_grped.first()
# which input compounds have cid match
has_cid = newPCL['compound name'].isin(mtch_first.index)

# reorder cid matches 
mtch_first = mtch_first.reindex(newPCL[has_cid]['compound name'])
newPCL['first_CID'] = np.nan
newPCL.ix[has_cid,'first_CID'] = mtch_first.values # add to dataframe

# write drug-cid pairs to a file
CID_ser = newPCL[has_cid]['first_CID']
CID_ser.index = newPCL[has_cid]['compound name']
CID_ser.name = 'pubchem_CID'
# oFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs_cids_list.csv'
oFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_cids_only.csv'
CID_ser.to_csv(oFile,sep=',',header=False,index=False,float_format='%.0f')

#########################
### add SMILE strings ### 
#########################

# load smile strings
sFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_SMILE.txt'
smileFrm = pd.read_csv(sFile,sep='\t',names=['cid', 'SMILE'])
smileSer = smileFrm.SMILE
smileSer.index = smileFrm.cid
smileSer = smileSer.reindex(newPCL.first_CID)
newPCL['SMILE'] = smileSer.values
# write new sheet to file
oFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs_SMILE.csv'
newPCL.to_csv(oFile,sep=',',index=False)

#########################
### add inCHI strings ### 
#########################

#load older most recent version of doc
oFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs_SMILE.csv'
pclFrm = pd.read_csv(oFile,sep=',')
#load and combine inchi keys from pubchem search
sFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/pubchem_inchi_results.txt'
inchFrm = pd.read_csv(sFile,sep='\t',names=['cid', 'InChI'])
inchSer = inchFrm.InChI
inchSer.index = inchFrm.cid
inchSer = inchSer.reindex(pclFrm.first_CID)
pclFrm['InChI'] = inchSer.values
# write new sheet to file
# oFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs_InChi.csv'
# pclFrm.to_csv(oFile,sep=',',index=False)


#########################
### inCHI to BRD number ### 
#########################

pFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/pert_id_inchikey.csv'
brdFrm = pd.read_csv(pFile)





