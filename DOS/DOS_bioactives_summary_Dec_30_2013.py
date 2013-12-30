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
import cmap.analytics.summly_null as SN
from statsmodels.distributions import ECDF

wkdir = '/xchip/cogs/projects/DOS/bioactivity_summary_Dec182013'
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
def load_summly_independent(iGold,mtrxSummly,index_row_by_pert_type=False):
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
    if index_row_by_pert_type:
        rowPertType = IST.get_row_meta('pert_type')
        rowPertIDs = IST.get_row_meta('id')
        iZip = zip(*[rowPertType,rowPertIDs])
        mRow = pd.MultiIndex.from_tuples(iZip, names=['pert_type','pert_id'])
        inSum.index = mRow
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

def median_conn_by_pert_type(inSum,dmsoSum,matrixType,rnkpt_thresh=90,graph=True):
    "1) what is the number of connections above a given rnkpt threshold \
    -what is the median for each unique pert"
    for pType in set(inSum.index.get_level_values('pert_type')): #loop through each pert_type
        pSum = inSum.ix[pType] #set matrix to summly results for a single pert_type
        pDMSO = dmsoSum.ix[pType]
        passMask = np.zeros_like(pSum.values)
        passMask[np.where(pSum.values > rnkpt_thresh)] = 1
        # count connections passed theshold
        passSum = np.sum(passMask,axis=0)
        passSer = pd.Series(data=passSum,index=pSum.columns)
        passSer.name = 'number_of_connections_pass_' + str(rnkpt_thresh) + '_rnkpt'
        passGrped = passSer.groupby(level='pert_id')
        dosMedConnect = passGrped.median()
        dosMedConnect.name = 'median_number_of_connections_above_' + str(rnkpt_thresh) + '_rnkpt'
        # repeat calculation for DMSOs
        passMaskDMSO = np.zeros_like(pDMSO.values)
        passMaskDMSO[np.where(pDMSO.values > rnkpt_thresh)] = 1
        passSumDMSO = np.sum(passMaskDMSO,axis=0)
        dmsoSer = pd.Series(data=passSumDMSO,index=pDMSO.columns)
        dmsoSer.name = 'number_of_connections_above_' + str(rnkpt_thresh) + '_rnkpt'
        if graph:
            min1 = np.min([np.min(dmsoSer.values),np.min(dosMedConnect)])
            max1 = np.max([np.max(dmsoSer.values),np.max(dosMedConnect)])
            h1 = plt.hist(dmsoSer,30,color='b',range=[min1,max1],label=['DMSO n=' + str(len(dmsoSer))],alpha=.4,normed=True)
            h3 = plt.hist(dosMedConnect,30,color='r',range=[min1,max1],label=['DOS n=' + str(len(dosMedConnect))],alpha=.3,normed=True) #
            plt.legend()
            plt.ylabel('normed freq',fontweight='bold')
            plt.xlabel('median counts ('+ matrixType + ' > ' + str(rnkpt_thresh) + ')',fontweight='bold')
            plt.title(pType + ' median connections (compounds collapsed by pert_id)')
            outF = os.path.join(wkdir, pType + '_median_summly_counts_pass_threshold.png')
            plt.savefig(outF, bbox_inches='tight',dpi=200)
            plt.close()
        # return dosMedConnect, dmsoSer

###################################
## summly connection consistency ##
###################################

def connection_overlap_median(inSum,dmsoSum,matrixType,nTop_connections=50,graph=True):
    "how consistent are the observed connections across cell lines? \
    -For all signatures of a given compound, calculate the pairwise \
    overlap among summly results \
    -take the median for each compound \
    -compare to DMSO overlap results "
    def get_top_set(x,nTop=nTop_connections):
        'make a set of the top summly connections'
        iTop = x[x<nTop]
        iList = set(iTop.index.values)
        return iList
    for pType in set(inSum.index.get_level_values('pert_type')): #loop through each pert_type
        pSum = inSum.ix[pType] #set matrix to summly results for a single pert_type
        pDMSO = dmsoSum.ix[pType]
        ### calculate median overlap for compound data
        rankedFrm = pSum.rank(axis=0,ascending=False)
        rankedFrm.columns = rankedFrm.columns.get_level_values('sig_id')
        topConnections = rankedFrm.apply(get_top_set,axis=0)
        if (topConnections.index == pSum.columns.get_level_values('sig_id')).all():
            topConnections.index = pSum.columns
        tcGrped = topConnections.groupby(level='pert_id')
        matrixSer = tcGrped.apply(test_overlap)
        overlapMed = matrixSer.apply(np.median)
        overlapMed = overlapMed[~np.isnan(overlapMed)]
        ### calculate median overlap for DMSO data
        #make DMSO groupings of the same size
        idGrped = pSum.groupby(level='pert_id',axis=1)
        pairSize = idGrped.size() # signature counts for each compound
        pairSize = pairSize[pairSize > 1]
        tmpFrm = pSum.reindex(columns=pairSize.index,level='pert_id')
        # make random dmso pairings to match compound pairing
        dmsoSigs = np.random.choice(pDMSO.columns.values,size=tmpFrm.shape[1])
        iZip = zip(*[tmpFrm.columns.get_level_values('pert_id'),dmsoSigs])
        mRow = pd.MultiIndex.from_tuples(iZip, names=['dos_pert_id','dmso_sig_id'])
        # rank DMSO summly results and calcualte overlap 
        rankedDmso = pDMSO.rank(axis=0,ascending=False)
        topDMSOConn = rankedDmso.apply(get_top_set,axis=0)
        mtchDMSOtop = topDMSOConn.reindex(dmsoSigs)
        mtchDMSOtop.index = mRow # match DMSO index to be like compounds
        tcDMSO = mtchDMSOtop.groupby(level='dos_pert_id')
        dmSer = tcDMSO.apply(test_overlap)
        oMedDMSO = dmSer.apply(np.median)
        oMedDMSO = oMedDMSO[~np.isnan(oMedDMSO)]
        if graph:
            min1 = np.min([np.min(oMedDMSO),np.min(overlapMed)])
            max1 = np.max([np.max(oMedDMSO),np.max(overlapMed)])
            h1 = plt.hist(oMedDMSO,30,color='b',range=[min1,max1],label=['DMSO'],alpha=.4,normed=True)
            h3 = plt.hist(overlapMed,30,color='r',range=[min1,max1],label=['DOS'],alpha=.3,normed=True) #
            plt.legend()
            plt.xlabel(pType + ' - median summly connection overlap (out of 50)')
            plt.ylabel('freq')
            plt.title('connection consistency across signatures')
            outF = os.path.join(wkdir, pType +  '_median_summly_connection_consistency.png')
            plt.savefig(outF, bbox_inches=0)
            plt.close()
def test_overlap(xSer):
    'test the set overlap among items in a Series \
    -return set of all pairwise overlap values'
    ns = len(xSer)
    if ns >1:
        overlapMtrx = np.zeros([ns,ns])
        for i1, ix1 in enumerate(xSer.index):
            s1 = xSer[ix1]
            if isinstance(s1,pd.Series):
                s1 = s1[0]
            for i2, ix2 in enumerate(xSer.index):
                s2 = xSer[ix2]
                if isinstance(s2,pd.Series):
                    s2 = s2[0]
                nO = len(s1.intersection(s2))
                overlapMtrx[i1,i2] = nO
                overlapMtrx[i2,i1] = nO                
        iUp = np.tril_indices(ns)
        overlapMtrx[iUp] = np.nan
        overlapArray = overlapMtrx.flatten()
        overlapArray = overlapArray[~np.isnan(overlapArray)]
        return overlapArray
    else:
        return np.zeros(0)

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

##########################################
### run functions to retrieve DOS info ###
##########################################


#get dmso results
sn = SN.SummNull(out=wkdir)
sn.load_dmso_summ_results(index_row_by_pert_type=True)
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
inSum,outSum = load_summly_independent(iGold,mtrxSummly,index_row_by_pert_type=True)
#get median rnkpt of n_top connections
# medianSer = median_n_connections(inSum,n_top=25)
# passSer = summ_connections_pass_thresh(inSum,outSum,rnkpt_thresh=90,graph=True)
# passSer = summ_conn_pass_thresh_vs_DMSO(inSum,sn.dmsoFrm,rnkpt_thresh=90,graph=True)
# dosMedConnect, dmsoSer = median_summ_conn_pass_thresh(inSum,
#                                                 outSum,
#                                                 sn.dmsoFrm,
#                                                 matrixType,
#                                                 rnkpt_thresh=90,
#                                                 graph=True)
# median_conn_by_pert_type(inSum,sn.dmsoFrm,matrixType,rnkpt_thresh=90,graph=True)
# resPreCalc = '/xchip/cogs/projects/connectivity/introspect/introspect_connectivity.txt' #
# specFrm, dosFrm = dos_introspect(resPreCalc,graph_metric='median_rankpt',graph=True)
connection_overlap_median(inSum,sn.dmsoFrm,matrixType,nTop_connections=50,graph=True)
# #save results to file
# outF = os.path.join(wkdir, 'DOS_signatures_counts_above_90_mrp4.txt')
# passSer.to_csv(outF,index=True,header=True,sep='\t')
# outF = os.path.join(wkdir, 'DOS_compounds_median_counts_above_90_mrp4.txt')
# dosMedConnect.to_csv(outF,index=True,header=True,sep='\t')
# outF = os.path.join(wkdir, 'DMSO_signatures_counts_above_90_mrp4.txt')
# dmsoSer.to_csv(outF,index=True,header=True,sep='\t')

# #seperate out by pert_type: 1) trt_cp vs. trt_sh
# #seperate out by pert_id, look at only unique compounds



## what is the relationship between median overlap and median CC or SS

