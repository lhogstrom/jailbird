'''
-identify the bioactive DOS compounds in CMAP
-calculate summary stats
-prep figure for manuscript

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
import cmap.analytics.summly_null as SN
from statsmodels.distributions import ECDF

wkdir = '/xchip/cogs/hogstrom/scripts/communique/Manus/M1/figure_code/figure_8_4'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

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
    # read in non-dos results
    ISO = gct.GCT(mtrxSummly)
    ISO.read(col_inds=list(iNonDos))
    outSum = ISO.frame
    gctSigs = ISO.get_column_meta('sig_id')
    gctPertIDs = ISO.get_column_meta('pert_id')
    # outSum.columns = gctSigs
    iZip = zip(*[gctPertIDs,gctSigs])
    mCol = pd.MultiIndex.from_tuples(iZip, names=['pert_id','sig_id'])
    outSum.columns = mCol
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

def summ_connections_pass_thresh(inSum,outSum,rnkpt_thresh=90,graph=True):
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
        plt.xlabel('counts ('+ matrixType + ' > ' + str(rnkpt_thresh) + ')',fontweight='bold')
        plt.title('DOS compounds with summly connections pass rnkpt ' + str(rnkpt_thresh))
        outF = os.path.join(wkdir, 'DOS_nonDOS_summly_counts_pass_threshold.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        plt.close()
    return passSer

def summ_conn_pass_thresh_vs_DMSO(inSum,dmsoSum,rnkpt_thresh=90,graph=True):
    "1) what is the number of connections past a given rnkpt threshold \
    how do DOS compounds compare to DMSO in summly space?"
    passMask = np.zeros_like(inSum.values)
    passMask[np.where(inSum.values > rnkpt_thresh)] = 1
    # count connections passed theshold
    passSum = np.sum(passMask,axis=0)
    passSer = pd.Series(data=passSum,index=inSum.columns)
    passSer.name = 'number_of_connections_pass_' + str(rnkpt_thresh) + '_rnkpt'
    passSer.index.name = 'sig_id'
    # repeat calculation for DMSOs
    passMaskNon = np.zeros_like(dmsoSum.values)
    passMaskNon[np.where(dmsoSum.values > rnkpt_thresh)] = 1
    passSumNon = np.sum(passMaskNon,axis=0)
    if graph:
        min1 = np.min([np.min(passSer.values),np.min(passSumNon)])
        max1 = np.max([np.max(passSer.values),np.max(passSumNon)])
        h1 = plt.hist(passSumNon,30,color='b',range=[min1,max1],label=['DMSO n=' + str(len(passSumNon))],alpha=.4,normed=True)
        h2 = plt.hist(passSer,30,color='r',range=[min1,max1],label=['DOS n=' + str(passSer.shape[0])],alpha=.3,normed=True) #
        plt.legend()
        plt.ylabel('normed freq',fontweight='bold')
        plt.xlabel('counts ('+ matrixType + ' > ' + str(rnkpt_thresh) + ')',fontweight='bold')
        plt.title('DOS signatures with summly connections pass rnkpt ' + str(rnkpt_thresh))
        outF = os.path.join(wkdir, 'DOS_DMSO_summly_counts_pass_threshold.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        plt.close()
    return passSer

def median_summ_conn_pass_thresh(inSum,outSum,dmsoSum,matrixType,rnkpt_thresh=90,graph=True):
    "1) what is the number of connections past a given rnkpt threshold \
    -what is the median for each unique pert"
    passMask = np.zeros_like(inSum.values)
    passMask[np.where(inSum.values > rnkpt_thresh)] = 1
    # count connections passed theshold
    passSum = np.sum(passMask,axis=0)
    passSer = pd.Series(data=passSum,index=inSum.columns)
    passSer.name = 'number_of_connections_pass_' + str(rnkpt_thresh) + '_rnkpt'
    passGrped = passSer.groupby(level='pert_id')
    dosMedConnect = passGrped.median()
    dosMedConnect.name = 'median_number_of_connections_above_' + str(rnkpt_thresh) + '_rnkpt'
    # repeat calculation for DMSOs
    passMaskDMSO = np.zeros_like(dmsoSum.values)
    passMaskDMSO[np.where(dmsoSum.values > rnkpt_thresh)] = 1
    passSumDMSO = np.sum(passMaskDMSO,axis=0)
    dmsoSer = pd.Series(data=passSumDMSO,index=dmsoSum.columns)
    dmsoSer.name = 'number_of_connections_above_' + str(rnkpt_thresh) + '_rnkpt'
    # repeat calculation for non-dos compounds
    passMaskNon = np.zeros_like(outSum.values)
    passMaskNon[np.where(outSum.values > rnkpt_thresh)] = 1
    passSumNon = np.sum(passMaskNon,axis=0)
    nonSer = pd.Series(data=passSumNon,index=outSum.columns)
    nonSer.name = 'number_of_connections_pass_' + str(rnkpt_thresh) + '_rnkpt'
    nonSer.index.name = 'sig_id'
    nonGrped = nonSer.groupby(level='pert_id')
    nonMedConnect = nonGrped.median()    
    if graph:
        min1 = np.min([np.min(passSer.values),np.min(passSumNon),np.min(dmsoSer.values)])
        max1 = np.max([np.max(passSer.values),np.max(passSumNon),np.min(dmsoSer.values)])
        h1 = plt.hist(dmsoSer,30,color='b',range=[min1,max1],label=['DMSO n=' + str(len(dmsoSer))],alpha=.4,normed=True)
        # h2 = plt.hist(nonMedConnect,30,color='g',range=[min1,max1],label=['non_DOS n=' + str(len(nonMedConnect))],alpha=.4,normed=True)
        h3 = plt.hist(dosMedConnect,30,color='r',range=[min1,max1],label=['DOS n=' + str(len(dosMedConnect))],alpha=.3,normed=True) #
        plt.legend()
        plt.ylabel('normed freq',fontweight='bold')
        plt.xlabel('median counts ('+ matrixType + ' > ' + str(rnkpt_thresh) + ')',fontweight='bold')
        plt.title('median connections (compounds collapsed by pert_id) - pass rnkpt ' + str(rnkpt_thresh))
        outF = os.path.join(wkdir, 'median_summly_counts_pass_threshold.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        plt.close()
        ### make cdf graph ####
        vals = np.linspace(min1,max1,100)
        dosEcdf = ECDF(dosMedConnect)
        dmsoEcdf = ECDF(dmsoSer)
        nonEcdf = ECDF(dosMedConnect)
        obsDos = dosEcdf(vals)
        obsDmso = dmsoEcdf(vals)
        obsNon = nonEcdf(vals)
        a1 = plt.plot(vals,obsDos,color='b',label=['DOS n=' + str(len(dosMedConnect))])
        a2 = plt.plot(vals,obsNon,color='g',label=['non_DOS n=' + str(len(nonMedConnect))])
        a3 = plt.plot(vals,obsDmso,color='r',label=['DMSO n=' + str(len(dmsoSer))]) #
        # plt.legend()
        plt.ylabel('F(x)',fontweight='bold')
        plt.xlabel('median counts ('+ matrixType + ' > ' + str(rnkpt_thresh) + ')',fontweight='bold')
        # plt.title('median connections pass rnkpt ' + str(rnkpt_thresh))
        outF = os.path.join(wkdir, 'median_summly_counts_cdf.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        plt.close()
    return dosMedConnect, dmsoSer

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
#get dmso results
sn = SN.SummNull(out=wkdir)
sn.load_dmso_summ_results()
dosBrds = get_dos_BRDs()
dosQuery, countSer = get_dos_signatures(dosBrds)
dosGold, countSerGold = get_dos_gold_signatures(dosBrds)
### use lass matrix
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_lass_n39560x7147.gctx'
# matrixType = 'rnkpt_indp_lass'
### use mrp4 mtrx
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_mrp4_n39560x7147.gctx'
matrixType = 'mrp4'
iGold = get_summly_dos_indeces(dosGold,mtrxSummly)
inSum,outSum = load_summly_independent(iGold,mtrxSummly)
#get median rnkpt of n_top connections
medianSer = median_n_connections(inSum,n_top=25)
# passSer = summ_connections_pass_thresh(inSum,outSum,rnkpt_thresh=90,graph=True)
passSer = summ_conn_pass_thresh_vs_DMSO(inSum,sn.dmsoFrm,rnkpt_thresh=90,graph=True)
dosMedConnect, dmsoSer = median_summ_conn_pass_thresh(inSum,
                                                outSum,
                                                sn.dmsoFrm,
                                                matrixType,
                                                rnkpt_thresh=90,
                                                graph=True)
# resPreCalc = '/xchip/cogs/projects/connectivity/introspect/introspect_connectivity.txt' #
# specFrm, dosFrm = dos_introspect(resPreCalc,graph_metric='median_rankpt',graph=True)

#save results to file
outF = os.path.join(wkdir, 'DOS_signatures_counts_above_90_mrp4.txt')
passSer.to_csv(outF,index=True,header=True,sep='\t')
outF = os.path.join(wkdir, 'DOS_compounds_median_counts_above_90_mrp4.txt')
dosMedConnect.to_csv(outF,index=True,header=True,sep='\t')
outF = os.path.join(wkdir, 'DMSO_signatures_counts_above_90_mrp4.txt')
dmsoSer.to_csv(outF,index=True,header=True,sep='\t')

#seperate out by pert_type: 1) trt_cp vs. trt_sh
#seperate out by pert_id, look at only unique compounds
