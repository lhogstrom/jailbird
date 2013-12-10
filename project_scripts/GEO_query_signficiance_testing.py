import getpass
import pymongo
import numpy as np, pandas as pd
from cmap.util import debug
import cmap
from os import path
import pandas as pd
import codecs
import cmap.io.gct as gct
import matplotlib.pyplot as plt

from pymongo import MongoClient
from cmap.util.mongo_utils import CredentialsObject
import cmap.util.mongo_utils as mu
from statsmodels.stats.multitest import fdrcorrection as fdr

wkdir = '/xchip/cogs/projects/GEO_summly'

### load GEO hypotheses being tested ####
CO = CredentialsObject()
CO.get_credentials()
server='23.23.153.188:27017,107.22.174.155:27017'
c = MongoClient(host=server) #,document_class='dex'
db = c['dex']
db.authenticate(CO.user,CO.pw)
# db.database['sig_info']
# collection = db['targets']
collection = db['summly_tags']
# collection.find_one()

g = collection.find()
# g = collection.find({'gene': 'HDAC1'})
dictList = list(g)
rFrame = pd.DataFrame(dictList)
rFrame.index = rFrame['_id']

#load in summly results
summFile = '/xchip/cogs/projects/GEO_summly/ext_summlies_mod.gctx'
gt = gct.GCT()
gt.read(summFile)
sFrm = gt.frame
#load in annotations file
annotFile = '/xchip/cogs/projects/GEO_summly/summlies.txt'
aFrm = pd.io.parsers.read_csv(annotFile,sep='\t')
aFrm.index = aFrm['_id']
#get expected connections for a geo query
summIDs = list(aFrm['_id'].values)
# rFrame[rFrame.source == aFrm['_id'][3]]
rFrame.reindex(index=summIDs)
matchFrm = rFrame[rFrame.source.isin(summIDs)]

# get pert_type for summly space sig_ids
summSpace = list(sFrm.index)
mc = mu.MongoContainer()
pIDfrm = mc.pert_info.find({'pert_id':{'$in':summSpace}},toDataFrame=True)

cpPerts = pIDfrm[pIDfrm.pert_type == 'trt_cp']['pert_id']
oePerts = pIDfrm[pIDfrm.pert_type == 'trt_oe']['pert_id']
shPerts = pIDfrm[pIDfrm.pert_type == 'trt_sh.cgs']['pert_id']
cpF = sFrm.reindex(index=cpPerts)
oeF = sFrm.reindex(index=oePerts)
shF = sFrm.reindex(index=shPerts)

#loop through the summly result for each GEO query
#calculate p-value of expected connections
geoNoSummly = []
geoNoRes = []
resDict = {}
for igeoID in summIDs:
    geoID = aFrm.ix[igeoID,'series']
    aName = aFrm.ix[igeoID,'a_name']
    if igeoID in rFrame.source.values:
        qFrm = rFrame[rFrame.source == igeoID]
        geoList.append(igeoID)
    else:
        geoNoRes.append(igeoID)
        continue        
    pertIDs = qFrm['sum_id'].values
    lmID = igeoID + '_lm_summly.txt'
    if lmID in sFrm.columns:
        cpRank, cpPercent, oeRank, oePercent, shRank, shPercent = calc_ranks(cpF,oeF,shF,lmID)
    else:
        geoNoSummly.append(lmID)
        continue
    for ipID in qFrm.index:
        eSer = qFrm.ix[ipID,:]
        pID = eSer['sum_id']
        pType = eSer['pert_type']
        eDir = eSer['direction']
        if pType == 'trt_cp':
            ePerc = cpPercent[pID]
            eRank = cpRank[pID]
            eRnkpt = cpF.ix[pID,lmID]
        if pType == 'trt_sh.cgs':
            ePerc = shPercent[pID]
            eRank = shRank[pID]
            eRnkpt = shF.ix[pID,lmID]
        if pType == 'trt_oe':
            ePerc = oePercent[pID]
            eRank = oeRank[pID]
            eRnkpt = oeF.ix[pID,lmID]
        #make table of expected connections
        #GEO, pert_id (expected), pert_type, rank, percent (p-value)
        resDict[geoID + ' - ' + pID] = {'geo_id':geoID,
                                'expected_connection':pID,
                                'pert_type':pType,
                                'rnkpt4':eRnkpt,
                                'rank':eRank,
                                'perc_rank_within_pert_type':ePerc,
                                'expected_direc':eDir,
                                'a_name':aName}
expectedFrm = pd.DataFrame(resDict)
expectedFrm = expectedFrm.T
expectedFrm.index.name = 'connection_pair'

expectedFrm['percent_rank_by_direction'] = np.nan
isPos = expectedFrm['expected_direc'] == 'Positively connected'
isNeg = expectedFrm['expected_direc'] == 'Negatively connected'
expectedFrm.ix[isPos,'percent_rank_by_direction'] = expectedFrm['perc_rank_within_pert_type']
expectedFrm.ix[~isPos,'percent_rank_by_direction'] = 1-expectedFrm['perc_rank_within_pert_type']

# plt.hist(expectedFrm['perc_rank'],30)
#check FDR using percent ranks
boolFDR, valFDR = fdr(pvals=expectedFrm['perc_rank_within_pert_type'].values) #,alpha=.05
expectedFrm['pass_FDR_with_percent_rank'] = boolFDR

### write expected connection summary
outF = wkdir+'/expected_connection_summary.txt'
expectedFrm.to_csv(outF,sep='\t',header=True,index=True)

## plot percent ranks overall
plt.hist(expectedFrm['perc_rank_within_pert_type'],30)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('percent rank of expected connections',fontweight='bold')
plt.title('All expected GEO connections')
outF = path.join(wkdir, 'percent_rank_expected_connections.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

## plot percent ranks accounting for expected direction
plt.hist(expectedFrm['percent_rank_by_direction'],30)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('percent rank of expected connections',fontweight='bold')
plt.title('expected GEO connections - accounting for direction')
outF = path.join(wkdir, 'percent_rank_by_direction_expected_connections.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

## plot percent ranks CP
plt.hist(expectedFrm[expectedFrm['pert_type'] == 'trt_cp']['perc_rank_within_pert_type'],30)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('percent rank of expected connections',fontweight='bold')
plt.title('CP expected GEO connections')
outF = path.join(wkdir, 'compounds_percent_rank_expected_connections.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

## plot percent ranks SH
plt.hist(expectedFrm[expectedFrm['pert_type'] == 'trt_sh.cgs']['perc_rank_within_pert_type'],30)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('percent rank of expected connections',fontweight='bold')
plt.title('CGS KD expected GEO connections')
outF = path.join(wkdir, 'sh_percent_rank_expected_connections.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

## plot percent ranks OE
plt.hist(expectedFrm[expectedFrm['pert_type'] == 'trt_oe']['perc_rank_within_pert_type'],30)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('percent rank of expected connections',fontweight='bold')
plt.title('OE expected GEO connections')
outF = path.join(wkdir, 'OE_percent_rank_expected_connections.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

#write files that didn't have a expected connections or summly results
noSumm = pd.Series(geoNoSummly)
noSumm.name = 'no_summly_results'
outF = wkdir+'/no_summly_results.txt'
noSumm.to_csv(outF,sep='\t',header=True,index=False)
#
noRes = pd.Series(geoNoRes)
noRes.name = 'no_expected_connections'
outF = wkdir+'/no_expected_connections.txt'
noRes.to_csv(outF,sep='\t',header=True,index=False)

def calc_ranks(cpF,oeF,shF,lmID):
    'calculate the rank and percent rank of summly results'
    cpS = cpF[lmID] #summly results for that GEO ID
    cpS = cpS.replace(-666,np.nan)
    cpS = cpS[~np.isnan(cpS)]
    cpS = cpS.order(ascending=False) # what do I do about 
    cpRank = pd.Series(data=np.arange(cpS.shape[0])+1,index=cpS.index)
    cpPercent = cpRank/float(len(cpRank))
    #oe
    oeS = oeF[lmID] #summly results for that GEO ID
    oeS = oeS.replace(-666,np.nan)
    oeS = oeS[~np.isnan(oeS)]
    oeS = oeS.order(ascending=False) # what do I do about 
    oeRank = pd.Series(data=np.arange(oeS.shape[0])+1,index=oeS.index)
    oePercent = oeRank/float(len(oeRank))
    #kd
    shS = shF[lmID] #summly results for that GEO ID
    shS = shS.replace(-666,np.nan)
    shS = shS[~np.isnan(shS)]
    shS = shS.order(ascending=False) # what do I do about 
    shRank = pd.Series(data=np.arange(shS.shape[0])+1,index=shS.index)
    shPercent = shRank/float(len(shRank))
    return cpRank, cpPercent, oeRank, oePercent, shRank, shPercent

### Obtain DMSOs 




