'''
-identify the bioactive DOS compounds in CMAP
-calculate summary stats
-prep for further analysis

Larson Hogstrom, 12/2013
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import cmap
import os
from matplotlib.ticker import NullFormatter

wkdir = '/xchip/cogs/projects/DOS/bioactivity_summary_Dec2013'

########################
## get DOS collection ##
########################

def get_dos_BRDs():
    'return Series of all DOS compounds - use the DOS icollection'
    #get all cps from DOS collection
    mc = mu.MongoContainer()
    pertInfo = mc.pert_info.find({'pert_icollection':'DOS'},
                {},toDataFrame=True)
    #check that it doesn't have a known pert_iname
    inameSer = pertInfo['pert_iname']
    inameFrm = pd.DataFrame(inameSer)
    #which values do not start with BRD?
    notBRDiname = pertInfo[~inameSer.str.contains('BRD')]
    isBRDiname = pertInfo[inameSer.str.contains('BRD')]
    dosBrds = isBRDiname['pert_id']
    return dosBrds

########################
## all DOS signatures ##
########################

def get_dos_signatures(dosBrds):
    "1) return signature info for all DOS compounds \
    2) return list counts of number of cell lines tested"
    CM = mu.CMapMongo()
    dosQuery = CM.find({'pert_id':{'$in':list(dosBrds)},'pert_type':'trt_cp'}, #, 
            {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
            toDataFrame=True)
    dosQuery.index = dosQuery['sig_id']
    dosSetLen = len(set(dosQuery['pert_id']))
    dosGrped = dosQuery.groupby(['pert_id'])
    countDict = {}
    for grp in dosGrped:
        grpName = grp[0]
        cellSet = set(grp[1]['cell_id'])
        nCells = len(cellSet)
        countDict[grpName] = nCells
    countSer = pd.Series(countDict)
    countMax = max(countSer)
    return dosQuery, countSer

#####################
## is_gold DOS ######
#####################

def get_dos_gold_signatures(dosBrds):
    "1) return signature info for all Gold DOS compounds \
    2) return list counts of number of cell lines tested"
    CM = mu.CMapMongo()
    dosGold = CM.find({'pert_id':{'$in':list(dosBrds)},'pert_type':'trt_cp','is_gold':True}, #, 
            {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
            toDataFrame=True)
    dosGold.index = dosGold['sig_id']
    dosGoldLen = len(set(dosGold['pert_id']))
    dosGrped = dosGold.groupby(['pert_id'])
    countGoldDict = {}
    for grp in dosGrped:
        grpName = grp[0]
        cellSet = set(grp[1]['cell_id'])
        nCells = len(cellSet)
        countGoldDict[grpName] = nCells
    countSerGold = pd.Series(countGoldDict)
    countSerGold.name = 'n_cell_lines'
    return dosGold, countSerGold

#which compounds are in independent mode space?
def get_summly_dos_indeces(dosGold,mtrxSummly):
    "1) obtain indices for all DOS compounds in the pre computed \
    summly matrix"
    gt = gct.GCT()
    gt.read_gctx_col_meta(mtrxSummly)
    gt.read_gctx_row_meta(mtrxSummly)
    indSummSigs = gt.get_column_meta('sig_id')
    indSummInames = gt.get_column_meta('pert_iname')
    sigSer = pd.Series(index=indSummSigs, data=indSummInames)
    #which dos cps are in the summly indpend
    summBrds = set(sigSer.index)
    goldDosBrds = set(dosGold['sig_id'].values)
    summGold = summBrds.intersection(goldDosBrds)
    # for the what is the median of the top n connections in summly independent mode?
    #get indices of gold DOS
    indSummSer = pd.Series(indSummSigs)
    indSer = pd.Series(index=indSummSer.values,data=indSummSer.index)
    iGold = indSer.reindex(list(summGold))
    return iGold

#load summly independent mode results
def load_summly_independent(iGold,mtrxSummly):
    "load dos compounds that have independent mode results - return dataframe"
    IST = gct.GCT(mtrxSummly)
    IST.read(col_inds=list(iGold.values))
    inSum = IST.frame
    gctSigs = IST.get_column_meta('sig_id')
    gctPertIDs = IST.get_column_meta('pert_id')
    # inSum.columns = gctSigs #index is just sig_id
    # hierarchical index - sig_id and pert_id
    iZip = zip(*[gctPertIDs,gctSigs])
    mCol = pd.MultiIndex.from_tuples(iZip, names=['pert_id','sig_id'])
    inSum.columns = mCol
    #read all non-dos summlies
    gt = gct.GCT()
    gt.read_gctx_col_meta(mtrxSummly)
    gt.read_gctx_row_meta(mtrxSummly)
    indSummSigs = gt.get_column_meta('sig_id')
    iNonDos = np.arange(len(indSummSigs))
    iNonDos = np.delete(iNonDos,iGold.values)
    # read in 
    ISO = gct.GCT(mtrxSummly)
    ISO.read(col_inds=list(iNonDos))
    outSum = ISO.frame
    gctSigs = ISO.get_column_meta('sig_id')
    outSum.columns = gctSigs
    return inSum, outSum #return dataframe of rankpt values

def median_n_connections(inSum,n_top=25):
    "what is the median rnkpt value of the n_top connections for each"
    resMtx = inSum.values.copy()
    resMtx.sort(axis=0)
    topVals = resMtx[-n_top:,:]
    topMedians = np.median(topVals,axis=0)
    medianSer = pd.Series(data=topMedians,index=inSum.columns)
    medianSer.name = 'median_of_top_25_connections'
    medianSer.index.name = 'sig_id'
    return medianSer

def summ_connections_pass_thresh(inSum,rnkpt_thresh=90,graph=True):
    "1) what is the number of connections past a given rnkpt threshold \
    how do DOS compounds compare to other compounds in summly space?"
    passMask = np.zeros_like(inSum.values)
    passMask[np.where(inSum.values > rnkpt_thresh)] = 1
    # count connections passed theshold
    passSum = np.sum(passMask,axis=0)
    passSer = pd.Series(data=passSum,index=inSum.columns)
    passSer.name = 'number_of_connections_pass_' + str(rnkpt_thresh) + '_rnkpt'
    passSer.index.name = 'sig_id'
    # repeat calculation for non-dos compounds
    passMaskNon = np.zeros_like(outSum.values)
    passMaskNon[np.where(outSum.values > rnkpt_thresh)] = 1
    passSumNon = np.sum(passMaskNon,axis=0)
    if graph:
        min1 = np.min([np.min(passSer.values),np.min(passSumNon)])
        max1 = np.max([np.max(passSer.values),np.max(passSumNon)])
        h1 = plt.hist(passSumNon,30,color='b',range=[min1,max1],label=['non_DOS n=' + str(len(passSumNon))],alpha=.4,normed=True)
        h2 = plt.hist(passSer,30,color='r',range=[min1,max1],label=['DOS n=' + str(passSer.shape[0])],alpha=.3,normed=True) #
        plt.legend()
        plt.ylabel('normed freq',fontweight='bold')
        plt.xlabel('counts (rnkpt_indp_lass > ' + str(rnkpt_thresh) + ')',fontweight='bold')
        plt.title('DOS compounds with summly connections pass rnkpt ' + str(rnkpt_thresh))
        outF = os.path.join(wkdir, 'DOS_summly_counts_pass_threshold.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        plt.close()
    return passSer

########################
## DOS sig_introspect ##
########################

### load sig_introspect precalculated result file
def dos_introspect(resPreCalc,graph_metric='median_rankpt',graph=True):
    "1) load pre-calculated sig_introspect results "
    specFrm =pd.io.parsers.read_csv(resPreCalc,sep='\t')
    specFrm.index = specFrm['pert_id']
    #which of these are dos cps?
    goldSet = countSerGold.index.values
    iSpecSet = set(specFrm['pert_id'])
    spectGold = iSpecSet.intersection(goldSet)
    #dos compounds for which we have introspect results
    dosFrm = specFrm[specFrm.pert_id.isin(spectGold)]
    if graph:
        iRnkpt = dosFrm[graph_metric]
        h1 = plt.hist(specFrm[graph_metric],30,color='b',range=[-80,100],label=['all_introspect_results'],alpha=.4,normed=True)
        h2 = plt.hist(iRnkpt,30,color='r',range=[-80,100],label='DOS_results',alpha=.3,normed=True)
        plt.legend()
        plt.ylabel('freq',fontweight='bold')
        plt.xlabel(graph_metric,fontweight='bold')
        # plt.title('summlySpace DOS compounds (' + str(overlapCount) + ') - cell lines is_gold')
        outF = os.path.join(wkdir, 'DOS_sig_introspect_'+ graph_metric+ '.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        plt.close()
    return specFrm, dosFrm # all results, dos only results

### run functions to retrieve DOS info 
dosBrds = get_dos_BRDs()
dosQuery, countSer = get_dos_signatures(dosBrds)
dosGold, countSerGold = get_dos_gold_signatures(dosBrds)
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_lass_n39560x7147.gctx'
iGold = get_summly_dos_indeces(dosGold,mtrxSummly)
inSum,outSum = load_summly_independent(iGold,mtrxSummly)
#get median rnkpt of n_top connections
medianSer = median_n_connections(inSum,n_top=25)
passSer = summ_connections_pass_thresh(inSum,rnkpt_thresh=90,graph=True)
resPreCalc = '/xchip/cogs/projects/connectivity/introspect/introspect_connectivity.txt' #
specFrm, dosFrm = dos_introspect(resPreCalc,graph_metric='median_rankpt',graph=True)


#how consistent are the observed connections across cell lines?
# overlap among top_n connections
n_top = 100
unBrd = set(inSum.columns.get_level_values(level=0)) #unique brds
for brd in unBrd:
    brdFrm = inSum[brd]
    nSigs = brdFrm.shape[1]
    if nSigs < 2: #skip if only one signature for that brd
        continue
    #what were to top n_connections for each signature
    brdFrm.idxmax()
    brdFrm

brdGrped = inSum.groupby(level='pert_id',axis=1)

def find_top(df, n_top=50, column='tip_pct'):
    # return df.sort_index(by=column)[-n:]
    ranked = df.rank(axis=0,ascending=False)
    mask = pd.DataFrame(data=np.zeros_like(df),index=df.index,columns=df.columns)
    mask[ranked<=n_top] = 1

    nTop = df[ranked<=n_top]

## which DOS compounds are in summly space?
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_sym_n7322x7322.gctx'
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_n7147x7147.gctx'
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_lass_n39560x7147.gctx'
gt = gct.GCT()
# gt.read_gctx_col_meta(mtrxSummly)
# gt.read_gctx_row_meta(mtrxSummly)
gt.read(mtrxSummly)
columnPerts = gt.get_column_meta('pert_id')
summFrm = gt.frame
summFrm.columns = columnPerts
# find dos cps in summly space
summBrds = summFrm.index.values
summSet = set(summBrds)
dosSet = set(countSer.index)
overlapSet = dosSet.intersection(summSet)
#plot hist of cell counts - dos in summly 
overlapSer = countSerGold.reindex(list(overlapSet))
overlapCount = len(overlapSer)
## plot
# countMax = max(overlapSer)
# bins = np.arange(countMax+1)
# plt.hist(overlapSer,bins)
# plt.ylabel('freq',fontweight='bold')
# plt.xlabel('number of cell lines',fontweight='bold')
# plt.title('summlySpace DOS compounds (' + str(overlapCount) + ') - cell lines is_gold')
# plt.xticks(bins)
# outF = os.path.join(wkdir, 'DOS_summly_cell_line_distribution.png')
# plt.savefig(outF, bbox_inches='tight',dpi=200)
# plt.close()



dosSer = sigSer.reindex(dosGold['sig_id'].values)


### make summary table:
# 1) pert_id
# 2) times_profiled_in_a2
# 3) times_gold_in_a2
# 4) is_gold_cell lines
goldGrped = dosGold.groupby('pert_id')
allGrped = dosQuery.groupby('pert_id')
#with vcap
sFrm = pd.DataFrame()
for brd in multiGold.index:
    allSigs = allGrped.groups[brd]
    nTested = len(allSigs)
    goldSigs = goldGrped.groups[brd]
    nGold = len(goldSigs)
    goldCells = dosGold.ix[goldSigs,'cell_id']
    goldCellSet = set(goldCells)
    noVCAP = goldCellSet.copy()
    if 'VCAP' in noVCAP:
        noVCAP.remove('VCAP')
    #which compounds are is gold in two cell lines OTHER than VCAP
    if len(noVCAP) > 1:
        nDict = {}
        nDict['times_profiled_in_a2'] = nTested
        nDict['times_gold_in_a2'] = nGold
        nDict['is_gold_cell_lines'] = list(goldCellSet)
        brdFrm = pd.DataFrame({brd:nDict})
        sFrm = pd.concat([sFrm,brdFrm],axis=1)
sFrm = sFrm.T
#convert list of cell lines to a comma seperated strings
sFrm['is_gold_cell_lines'] = sFrm['is_gold_cell_lines'].str.join(',')
sFrm = sFrm.reindex(columns=['times_profiled_in_a2','times_gold_in_a2','is_gold_cell_lines'])
outF = os.path.join(wkdir, 'DOS_gold_cell_line_summary.txt')
sFrm.to_csv(outF,sep='\t', header=True)

### for the 200 in summly space ---> run sig_introspect 
outIntrospect = wkdir + '/sig_introspect_DOS_summly_space_non_gold'
if not os.path.exists(outIntrospect):
    os.mkdir(outIntrospect)
#get all sig_ids for dos cps in summlyspace
CM = mu.CMapMongo()
summSpaceQuery = CM.find({'pert_id':{'$in':list(overlapSet)},'pert_type':'trt_cp'}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True,'distil_ss':True,'distil_cc_q75':True},
        toDataFrame=True)
qSer = summSpaceQuery['sig_id']
outF = outIntrospect + '/sig_ids_dos_summSpace.grp'
qSer.to_csv(outF,index=False,header=False)
#run sig_introspect
# cmd = ' '.join(['rum -q hour',
#      '-d sulfur_io=100',
cmd = ' '.join(['rum -q hour',
     '-o ' + outIntrospect,
     '-x sig_introspect_tool ',
     '--sig_id ' + outF,
     '--query_group pert_id',
     '--metric wtcs',
     '--out ' + outIntrospect])
os.system(cmd)
#which compounds have a pre-calculated introspect result?


# ### run summly for bioactives - first in matched mode for the 234, then in independent mode

## how do gold signatures correlate across cell lines?
multiRes = dosGold[dosGold.pert_id.isin(multiGold.index)]
# load Expression Signatures
sigs = multiRes['sig_id'].values
afPath = cmap.score_path
gt = gct.GCT()
gt.read(src=afPath,cid=list(sigs),rid='lm_epsilon')
expFrm = gt.frame

pairCorr = np.corrcoef(expFrm.values.T)
corrFrm = pd.DataFrame(data=pairCorr,index=expFrm.columns,columns=expFrm.columns)

multiRes.index = multiRes['sig_id']
pertGrped = multiRes.groupby('pert_id')
#coherence of signatures across all cell lines
corrDict = {}
randDict = {}
observedCorrs = []
randCorrs = []
for brd in pertGrped.groups:
    sigs = pertGrped.groups[brd]
    smFrm = corrFrm.reindex(index=sigs,columns=sigs)
    nm = smFrm.shape[0]
    iUp = np.tril_indices(nm)
    smFrm.values[iUp] = np.nan
    corrArray = smFrm.unstack().values
    corrArray = corrArray[~np.isnan(corrArray)]
    corrDict[brd] = corrArray
    observedCorrs.extend(corrArray)
    # null distribution --> random pairs of sig_ids
    iRand = np.random.choice(multiRes.shape[0],nm,replace=False)
    randSigs = multiRes.index[iRand]
    rndFrm = corrFrm.reindex(index=randSigs,columns=randSigs)
    rndFrm.values[iUp] = np.nan
    rndArray = rndFrm.unstack().values
    rndArray = rndArray[~np.isnan(rndArray)]
    randDict[brd] = rndArray
    randCorrs.extend(rndArray)
# make histogram of self-self connections vs. random 
h1 = plt.hist(observedCorrs,30,color='b',range=[-1,1],label=['observed'],alpha=.4)
h2 = plt.hist(randCorrs,30,color='r',range=[-1,1],label='random',alpha=.4)
plt.legend()
plt.xlabel('corr')
plt.ylabel('freq')
plt.title('DOS compound self connections across cell lines')
outF = os.path.join(wkdir, 'DOS_pairwise_self_connection_distribution.png')
plt.savefig(outF, bbox_inches=0)
plt.close()

FileWtcs = '/xchip/cogs/projects/connectivity/query/wtcs.lm50/sim_wtcs.lm50_COMBINED.gctx'
gt = gct.GCT()
gt.read(src=FileWtcs,cid=list(sigs),rid=list(sigs))
wtcsFrm = gt.frame


########################
## make SC Plots #######
########################

# x = dosGold['distil_cc_q75'].values
# y = dosGold['distil_ss'].values
# # y = nsample
# nullfmt   = NullFormatter()         # no labels
# # definitions for the axes
# left, width = 0.1, 0.65
# bottom, height = 0.1, 0.65
# bottom_h = left_h = left+width+0.02
# #
# rect_scatter = [left, bottom, width, height]
# rect_histx = [left, bottom_h, width, 0.2]
# rect_histy = [left_h, bottom, 0.2, height]
# # start with a rectangular Figure
# plt.figure(1, figsize=(8,8))
# axScatter = plt.axes(rect_scatter)
# axHistx = plt.axes(rect_histx)
# axHisty = plt.axes(rect_histy)
# # no labels
# axHistx.xaxis.set_major_formatter(nullfmt)
# axHisty.yaxis.set_major_formatter(nullfmt)
# # the scatter plot:
# axScatter.scatter(x, y, marker=".", alpha=.1) #,marker='.'
# binwidth = 0.25
# # xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
# xymax = np.max(np.fabs(y))
# # lim = ( int(xymax/binwidth) + 1) * binwidth
# lim = ( int(xymax/binwidth) + xymax) * binwidth
# #set 
# xMin = -.3
# xMax = 1
# xBinwidth = .025
# axScatter.set_xlim( (xMin, xMax) )
# axScatter.set_ylim( (0, lim) )
# #make hist
# binsX = np.arange(xMin, xMax + xBinwidth, xBinwidth)
# binsY = np.arange(0, lim + binwidth, binwidth)
# # axHistx.hist(x, bins=bins)
# axHistx.hist(x, bins=binsX, range=[xMin,xMax])
# axHisty.hist(y, bins=binsY, orientation='horizontal')
# axHistx.axis('off')
# axHisty.axis('off')
# #set lims
# axScatter.set_ylabel('ss',fontweight='bold')
# axScatter.set_xlabel('CC_q75',fontweight='bold')
# # axScatter.set_title('DMSO signatures - LINCS core cell lines (n = 6471)')
# outF = os.path.join(wkdir, 'DOS_SC_plot_all_is_gold.png')
# plt.savefig(outF, bbox_inches='tight',dpi=200)
# plt.close()

