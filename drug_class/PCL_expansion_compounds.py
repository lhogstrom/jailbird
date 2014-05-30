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
import cmap.io.gmt as gmt

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
### inCHI to InchiKey ### 
#########################

#load inchi key
iFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/CID-InChI-Key'
iKey = pd.read_csv(iFile,sep='\t',header=None,names=['CID','InChi-key'])

#load older most recent version of doc
oFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs_InChi.csv'
pclFrm = pd.read_csv(oFile,sep=',')
has_cid = ~np.isnan(pclFrm.first_CID)
cid_list = pclFrm.ix[has_cid,'first_CID']
iKey_match = iKey.reindex(cid_list.values)
pclFrm.ix[has_cid,'InChiKey'] = iKey_match['InChi-key'].values
oFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs_inChiKey.csv'
# pclFrm.to_csv(oFile,sep=',',index=False)

#############################
### InChiKey to BRD number ### 
#############################

#load older most recent version of doc
oFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs_inChiKey.csv'
pclFrm = pd.read_csv(oFile,sep=',')

pFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/pert_id_inchikey.csv'
brdFrm = pd.read_csv(pFile)

# check for inchikey and inchikey desalted
brdFrm = brdFrm[~brdFrm.INCHIKEY.isnull()] # drop rows without inchikey
key_match1 = brdFrm[brdFrm.INCHIKEY.isin(pclFrm.InChiKey)]
key_match2 = brdFrm[brdFrm.INCHIKEY_DESALTED.isin(pclFrm.InChiKey)]
iUnion = set.union(set(key_match1.index),set(key_match2.index))
key_match = brdFrm.reindex(iUnion)

# reindex to match summary sheet
BRDser = key_match.BROADID
BRDser.index = key_match.INCHIKEY
BRDser = BRDser.reindex(pclFrm.InChiKey)

### check if the BRDs are in CMAP ###
brdSet = set(BRDser.values)
mc = mu.MongoContainer()
cmapFrm = mc.pert_info.find({'pert_id':{'$in':list(brdSet)}},{},toDataFrame=True)

# add BRD IDs to summary sheet
pclFrm['Broad_ID'] = BRDser.values
pclFrm['cmap_brd_match'] = BRDser.isin(cmapFrm.pert_id).values

# pert_iname match to common name 
mc = mu.MongoContainer()
inameFrm = mc.pert_info.find({'pert_type':'trt_cp'},{'pert_iname':True,'pert_id':True},toDataFrame=True)
inameGrped = inameFrm.groupby('pert_iname')
inameSet = inameGrped.first()
inameSet['pert_iname'] = inameSet.index
inameLower = inameSet.copy()
inameLower['pert_iname'] = inameSet.pert_iname.str.lower()
#match to pcls
pclLower = pclFrm['compound name'].str.lower()
iname_match = pclLower.isin(inameLower['pert_iname'].values)
pclFrm['cmap_pert_iname_match'] = iname_match

### is_PCL (which compounds are already in a PCL)
cFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140221/cliques.gmt'
cliqueGMT = gmt.read(cFile)
cliqFrm = pd.DataFrame(cliqueGMT)
### list of compounds by cliques
brdCliq = pd.Series()
for icliq,cliq in enumerate(cliqFrm.desc):
    cliqMod = cliqFrm.ix[icliq,'desc']
    brds = cliqFrm.ix[icliq,'sig']
    for brd in brds:
        brdCliq[brd] = cliqMod

# which of these compounds already belong to a PCL?
pclBrds = pclFrm.Broad_ID[~pclFrm.Broad_ID.isnull()]
pclBrds = pclFrm.Broad_ID
is_pcl = pclBrds[pclBrds.isin(brdCliq.index)]
cliq_match = brdCliq.reindex(is_pcl.values)
cliq_match.index = is_pcl.index
pclFrm['is_PCL'] = np.nan
pclFrm.ix[cliq_match.index,'is_PCL'] = cliq_match.values

# write output to file
oFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/new_compounds_for_PCLs_BROADID.csv'
pclFrm.to_csv(oFile,sep=',',index=False)

