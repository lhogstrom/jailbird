'''
-use DMSO to access statstical significance of summly results
-examine false positive rates using ecdf of observed and dmso queries

Larson Hogstrom, 1/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import os
import cmap.analytics.summly_null as SN
from statsmodels.distributions import ECDF
import cmap.io.gmt as gmt
import cmap.util.progress as update
from matplotlib import cm
import multiprocessing
import time

wkdir = '/xchip/cogs/projects/connectivity/false_positive_rates/17Jan2014'
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

def get_summly_ind_compounds(dosGold,mtrxSummly):
    '''1) return non-DOS sig_ids of compounds that are in summly space \
    summly matrix

    Returns:
    pandas series - index = sig_ids, values = indices in summly matrix
    '''
    gt = gct.GCT()
    gt.read_gctx_col_meta(mtrxSummly)
    gt.read_gctx_row_meta(mtrxSummly)
    indSummSigs = gt.get_column_meta('sig_id')
    indSummPType = gt.get_column_meta('pert_type')
    indSummInames = gt.get_column_meta('pert_iname')
    # sigSer = pd.Series(index=indSummSigs, data=indSummInames)
    typeSer = pd.Series(index=indSummSigs, data=indSummPType)
    isCp = typeSer[typeSer == 'trt_cp']
    #which dos cps are in the summly indpend
    summBrds = set(isCp.index)
    goldDosBrds = set(dosGold['sig_id'].values)
    summGold = summBrds.difference(goldDosBrds)
    indSummSer = pd.Series(indSummSigs)
    indSer = pd.Series(index=indSummSer.values,data=indSummSer.index)
    iNonDos = indSer[indSer.index.isin(summGold)]
    # iNonDos = pd.Series(list(summGold))
    return iNonDos

def ecdf_calc(inSum,dmsoSum,matrixType,graph=True,fpr_max=True):
    '''
    -create empirical cdf for observed and dmso 

    Parameters:
    -----------

    '''
    #look at edcf by row
    seriesList = []
    progress_bar = update.DeterminateProgressBar('ecdf calculation')
    for ii,ix in enumerate(inSum.index):
        progress_bar.update('count', ii, len(inSum.index))
        pID = ix[1]
        obsVec = inSum.ix[ix]
        dmsoVec = dmsoSum.ix[ix]
        # flip sign of rnkpt values
        # evaluate ecdf
        oecdf = ECDF(obsVec)
        decdf = ECDF(dmsoVec)
        # min1 = np.min([np.min(obsVec),np.min(dmsoVec)])
        # max1 = np.max([np.max(obsVec),np.max(dmsoVec)])
        # vals = np.linspace(min1,max1,100)
        vals = np.linspace(-100,100,201)
        oEval = oecdf(vals)
        dEval = decdf(vals)
        # make individual plots
        # fdrVec = dEval / oEval
        fdrVec = (1 - dEval) / (1 - oEval) # looking for positive connections
        fdrSer = pd.Series(data=fdrVec, index=vals)
        if fpr_max:
            fdrSer[fdrSer >= 1] = 1
        fdrSer.name = ix
        seriesList.append(fdrSer)
        if graph:        
            fig = plt.figure(1, figsize=(10, 10))
            plt.subplot(2,1,1)
            a1 = plt.plot(vals,oEval,color='b',label='observed n=' + str(len(obsVec)))
            a3 = plt.plot(vals,dEval,color='r',label='DMSO n=' + str(len(dmsoVec))) #
            plt.legend(loc=2)
            plt.ylabel('F(x)',fontweight='bold')
            # plt.xlabel(matrixType,fontweight='bold')
            plt.title('ecdf for summly row - ' + pID)
            plt.subplot(2,1,2)
            h1 = plt.hist(obsVec,30,color='b',range=[-100,100],label=['observed'],alpha=.4,normed=True)
            h2 = plt.hist(dmsoVec,30,color='r',range=[-100,100],label='DMSO',alpha=.3,normed=True)
            # plt.legend()
            plt.ylabel('freq',fontweight='bold')
            plt.xlabel(matrixType,fontweight='bold')
            outF = os.path.join(wkdir, pID + '_ecdf.png')
            plt.savefig(outF, bbox_inches='tight',dpi=200)
            plt.close()
    fpFrame = pd.concat(seriesList, axis=1, keys=[s.name for s in seriesList])
    fpFrame = fpFrame.T
    mCol = pd.MultiIndex.from_tuples(fpFrame.index, names=['pert_type','pert_id'])
    fpFrame.index = mCol
    return fpFrame
    ### run this procedure without a loop using 'apply'
    # def ecdf_by_row(x):
    #     'make a set of the top summly connections'
    #     # iTop = x[x<nTop]
    #     # iList = set(iTop.index.values)
    #     # return iList
    # for pType in set(inSum.index.get_level_values('pert_type')): #loop through each pert_type
    #     fdrFrm = inSum.apply(get_top_set,axis=0)

def fdr_heatmaps(fdrFrm):
    '''
    given a DataFrame of FDR rates (given rnkpt threshold) - make heatmap

    '''
    # order acording to highest false positive rate @ rnkpt 90
    fpSort = fdrFrm.sort(90)
    # plot result
    fig = plt.figure(1, figsize=(10, 10))
    plt.imshow(fpSort.values,
        interpolation='nearest',
        aspect='auto',
        cmap=cm.gray_r)
    tickRange = range(0,fdrFrm.shape[1],20)
    xtcks = [str(x) for x in fpFrame.columns[tickRange]]
    plt.xticks(tickRange, xtcks)
    # plt.yticks(np.arange(len(ytcks)),ytcks)
    plt.colorbar()
    plt.xlabel(matrixType + ' threshold')
    plt.ylabel('unique perturbations')
    plt.title('summly false positive rate - based on DMSO')
    out = wkdir + '/false_positive_matrix_' + matrixType + '_threshold.png'
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    # heatmap by pert_type
    fpGrped = fdrFrm.groupby(level='pert_type')
    for grp in fpGrped.groups:
        grpFrm = fpGrped.get_group(grp)
        grpSort = grpFrm.sort(90)
        fig = plt.figure(1, figsize=(10, 10))
        plt.imshow(grpSort.values,
            interpolation='nearest',
            aspect='auto',
            cmap=cm.gray_r)
            # vmin=0, 
            # vmax=1,
        tickRange = range(0,fdrFrm.shape[1],20)
        xtcks = [str(x) for x in fpFrame.columns[tickRange]]
        plt.xticks(tickRange, xtcks)
        # plt.yticks(np.arange(len(ytcks)),ytcks)
        plt.colorbar()
        plt.xlabel(matrixType + ' threshold')
        plt.ylabel('unique perturbations')
        plt.title(grp +' summly false positive rate - based on DMSO')
        out = wkdir + '/' + grp + '_false_positive_matrix_' + matrixType + '_threshold.png'
        plt.savefig(out, bbox_inches='tight')
        plt.close()

#get dmso results
sn = SN.SummNull(out=wkdir)
sn.load_dmso_summ_results(index_row_by_pert_type=True,summly_type='mrp4')
# sn.load_dmso_summ_results(index_row_by_pert_type=True,summly_type='lass')
dosBrds = get_dos_BRDs()
dosGold, countSerGold = get_dos_gold_signatures(dosBrds)
### use lass matrix
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_lass_n39560x7147.gctx'
# matrixType = 'rnkpt_indp_lass'
### use mrp4 mtrx
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_mrp4_n39560x7147.gctx'
matrixType = 'mrp4'
iGold = get_summly_dos_indeces(dosGold,mtrxSummly)
# iNonDos = get_summly_ind_compounds(dosGold,mtrxSummly)
inSum,outSum = load_summly_independent(iGold,mtrxSummly,index_row_by_pert_type=True)
# fdr stats
fdrFrm = ecdf_calc(inSum,sn.dmsoFrm,matrixType,graph=False,fpr_max=True)
fdr_heatmaps(fdrFrm)

