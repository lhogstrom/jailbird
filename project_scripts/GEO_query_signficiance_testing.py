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
from matplotlib.ticker import NullFormatter

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
geoList = []
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
outF = wkdir+'/no_expected_connections_table.txt'
noExpectFrm = aFrm.reindex(geoNoRes)
noExpectFrm.to_csv(outF,sep='\t',header=True,index=False)

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

### Obtain DMSOs for permutation testing
CM = mu.CMapMongo()
cellList = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
#grab only well combined dmsos:
# dmsoQuery = CM.find({'pert_iname':'DMSO','cell_id':{'$in':cellList},'distil_nsample':{'$gte':2,'$lt':5}}, #, 
#         {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True,'brew_prefix':True,'up50_lm':True,'dn50_lm':True},
#         toDataFrame=True)
#grab all DMSOs:
dmsoQuery = CM.find({'pert_iname':'DMSO','cell_id':{'$in':cellList},'distil_cc_q75':{'$lt':.6},'distil_ss':{'$lt':7},'distil_nsample':{'$lt':7}}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True,'brew_prefix':True,'distil_nsample':True},
        toDataFrame=True)
# leave out signatures without distil_cc
is_null = dmsoQuery['distil_cc_q75'] == -666
dmsoQuery = dmsoQuery[~is_null.values]
# remove anything on from the CNS plates
is_cns = [x[:3] == 'CNS' for x in dmsoQuery['sig_id']]
dmsoQuery = dmsoQuery[~np.array(is_cns)]

# plateSet = set(dmsoQuery['brew_prefix'])
ccQ75 = dmsoQuery['distil_cc_q75']
SS = dmsoQuery['distil_ss']
# write sig_ids to file
outF = wkdir+'/dmso_sig_ids.grp'
dmsoQuery['sig_id'].to_csv(outF,sep='\t',header=False,index=False)

#plot relationship between CC and distil_nsample
nsample = dmsoQuery['distil_nsample']
plt.plot(nsample,ccQ75,'.')
plt.xlabel('distil_nsample',fontweight='bold')
plt.ylabel('CC_q75',fontweight='bold')
plt.title('DMSO signatures - LINCS core cell lines (n = 7182)')
outF = path.join(wkdir, 'LINCS_DMSOs_cc_by_distil_nsample.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()


# simple SC plot
# ccQ75 = dmsoQuery['distil_cc_q75']
# SS = dmsoQuery['distil_ss']
# plt.plot(ccQ75,SS,'.')
# plt.ylabel('ss',fontweight='bold')
# plt.xlabel('CC_q75',fontweight='bold')
# plt.title('DMSO signatures - LINCS core cell lines (n = 6471)')
# outF = path.join(wkdir, 'LINCS_DMSOs_SC_plot.png')
# plt.savefig(outF, bbox_inches='tight',dpi=200)
# plt.close()

### QQ plot with histogram
x = ccQ75
# y = SS
y = nsample
nullfmt   = NullFormatter()         # no labels
# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left+width+0.02
#
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]
# start with a rectangular Figure
plt.figure(1, figsize=(8,8))
axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)
# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)
# the scatter plot:
axScatter.scatter(x, y) #,marker='.'
# now determine nice limits by hand:
binwidth = 0.25
# xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
xymax = np.max(np.fabs(y))
# lim = ( int(xymax/binwidth) + 1) * binwidth
lim = ( int(xymax/binwidth) + xymax) * binwidth
#set 
xMin = -.3
xMax = 1
xBinwidth = .05
axScatter.set_xlim( (xMin, xMax) )
axScatter.set_ylim( (0, lim) )
#make hist
binsX = np.arange(xMin, xMax + xBinwidth, xBinwidth)
binsY = np.arange(0, lim + binwidth, binwidth)
# axHistx.hist(x, bins=bins)
axHistx.hist(x, 30, range=[xMin,xMax])
axHisty.hist(y, bins=binsY, orientation='horizontal')
#set lims
# axHistx.set_xlim( axScatter.get_xlim() )
# axHisty.set_ylim( axScatter.get_ylim() )
# plt.show()
# axScatter.set_ylabel('ss',fontweight='bold')
axScatter.set_ylabel('distil_nsample',fontweight='bold')
axScatter.set_xlabel('CC_q75',fontweight='bold')
# axScatter.set_title('DMSO signatures - LINCS core cell lines (n = 6471)')
# outF = path.join(wkdir, 'LINCS_DMSOs_SC_plot_hist.png')
outF = path.join(wkdir, 'LINCS_DMSOs_cc_by_distil_nsample.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

# make gmt for DMSO queries
gmtListUp = []
gmtListDn = []
for i1 in dmsoQuery.index:
    gmtDictUp = {}
    gmtDictUp['id'] = dmsoQuery.ix[i1,'sig_id']
    gmtDictUp['desc'] = dmsoQuery.ix[i1,'sig_id']
    gmtDictUp['sig'] = dmsoQuery.ix[i1,'up50_lm']
    gmtListUp.append(gmtDictUp)
    gmtDictDn = {}
    gmtDictDn['id'] = dmsoQuery.ix[i1,'sig_id']
    gmtDictDn['desc'] = dmsoQuery.ix[i1,'sig_id']
    gmtDictDn['sig'] = dmsoQuery.ix[i1,'dn50_lm']
    gmtListDn.append(gmtDictDn)
gmtOutUp = wkdir + '/DMSO_lincs_core_lines_Up_lm50.gmt'
gmtOutDn = wkdir + '/DMSO_lincs_core_lines_Dn_lm50.gmt'
gmt.write(gmtListUp,gmtOutUp)
gmt.write(gmtListDn,gmtOutDn)


larger_query_matrix = '/xchip/cogs/projects/connectivity/query/wtcs.lm50/sim_wtcs.lm50_COMBINED.gctx'
gt = gct.GCT()
gt.read(larger_query_matrix,)


#once the gmt is made:
# 1) run sig_query_tool
# 2) run sig_summly_tool
# 3) 



# metric = 'wtcs'
# queryDir = work_dir + '/ctd2_sig_query'
# if not os.path.exists(queryDir):
#     os.mkdir(queryDir)
# cmd = ' '.join(['rum -q local -f sig_query_tool',
#          '--sig_id ' + file1,
#          '--metric ' + metric,
#          '--column_space full',
#          '--out ' + queryDir,
#          '--mkdir false',
#          '--save_tail false'])
#          # '--row_space bing', 
# os.system(cmd)

# for cellLine in cellList:
#     # cellLine = 'A375'
#     wkdir = '/xchip/cogs/projects/NMF/BRAF_PCL/' + cellLine
#     if not os.path.exists(wkdir):
#         os.mkdir(wkdir)
#     # get signature annotations from cmap database
#     CM = mu.CMapMongo()
#     goldQuery = CM.find({'is_gold' : True,'pert_id':{'$in':brdAllGroups},'cell_id':cellLine,'pert_dose':{'$gt':1}}, #, 
#             {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
#             toDataFrame=True)
#     indRep = [x.replace(":",".") for x in goldQuery['sig_id']]
#     indRep = [x.replace("-",".") for x in indRep]
#     goldQuery.index = indRep
#     # goldQuery.index = goldQuery['sig_id']
#     dmsoQuery = CM.find({'pert_iname':'DMSO','cell_id':cellLine,'distil_nsample':{'$gte':2,'$lt':5}}, #, 
#             {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
#             toDataFrame=True)
#     indRep = [x.replace(":",".") for x in dmsoQuery['sig_id']]
#     indRep = [x.replace("-",".") for x in indRep]



# random dos cp lookup
# CM = mu.CMapMongo()
# cellList = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
# dosQuery = CM.find({'pert_id':'BRD-K77790153'}, #, 
#         {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True,'brew_prefix':True,'pert_dose':True},
#         toDataFrame=True)
# outF = wkdir+'/BRD-K77790153_signature_info.txt'
# dosQuery.to_csv(outF,sep='\t',header=True,index=False)



