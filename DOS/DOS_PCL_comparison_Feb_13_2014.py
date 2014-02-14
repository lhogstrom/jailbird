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
import os
import cmap.analytics.summly_null as SN
from statsmodels.distributions import ECDF
import cmap.io.gmt as gmt
import cmap.util.progress as update
from matplotlib import cm
import cmap.plot.colors as ccol

wkdir = '/xchip/cogs/projects/DOS/PCL_comparison_Feb162014'
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

############################
## summly matched DOS ######
############################

#which compounds are in independent mode space?
def get_summly_dos_indeces(dosBrds,mtrxSummly):
    "1) return dos compounds that are in summly matched space "
    gt = gct.GCT()
    gt.read(mtrxSummly)
    # indSummSigs = gt.get_column_meta('sig_id')
    # indSummInames = gt.get_column_meta('pert_iname')
    summFrm = gt.frame
    sigSer = pd.Series(index=summFrm.index, data=summFrm.columns)
    dosSer = sigSer[sigSer.index.isin(dosBrds)]
    return dosSer

##########################################
### run functions to retrieve DOS info ###
##########################################

dosBrds = get_dos_BRDs()
# dosQuery, countSer = get_dos_signatures(dosBrds)
# dosGold, countSerGold = get_dos_gold_signatures(dosBrds)
### use lass matrix
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_lass_n39560x7147.gctx'
# matrixType = 'rnkpt_indp_lass'
### use mrp4 mtrx
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_mrp4_n39560x7147.gctx'
# matrixType = 'mrp4'
### use lass matched matrix
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_n7147x7147.gctx'
matrixType = 'rnkpt_matched_lass'
dosSer = get_summly_dos_indeces(dosBrds,mtrxSummly)
outGRP = wkdir + '/dos_summly_ids.grp'
dosSer.to_csv(outGRP,index=False)

##########################################
### get brds for the whole matrix ###
##########################################

gt = gct.GCT()
gt.read(mtrxSummly)
summFrm = gt.frame
sigSer = pd.Series(index=summFrm.index, data=summFrm.columns)
outGRP = wkdir + '/summly_matched_ids.grp'
sigSer.to_csv(outGRP,index=False)

##########################################
### construct sig_cliquescore_tool command ###
##########################################

groupGMT = '/xchip/cogs/projects/pharm_class/rnwork/cliques/cpd_groups_n147.gmt'
cliqueGMT = gmt.read(groupGMT)
cliqFrm = pd.DataFrame(cliqueGMT)

cmd = ', '.join(['sig_cliquescore_tool(\'score\', \'' + mtrxSummly+ '\'',
         '\'summly_id\', \'' + outGRP + '\'',
         '\'clique\', \'' + groupGMT + '\'',
         '\'out\', \'' + wkdir + '\')'])

#########################################
### load sig_cliquescore_tool results ###
#########################################

# cliques against DOS compounds
cFile = '/xchip/cogs/projects/DOS/PCL_comparison_Feb162014/feb13/my_analysis.sig_cliquescore_tool.2014021313530491/clique_median_n145x234.gctx'
gt = gct.GCT()
gt.read(cFile)
cliqDos = gt.frame

# cliques against all compounds
cFileFull = '/xchip/cogs/projects/DOS/PCL_comparison_Feb162014/feb14/my_analysis.sig_cliquescore_tool.2014021411101091/clique_median_n145x7147.gctx'
gt2 = gct.GCT()
gt2.read(cFileFull)
cliqFull = gt2.frame
cliqFull = cliqFull.reindex(sigSer.values)
cliqFull.index = sigSer.index
# transpose matrix and write
# cT = cliqFull.T 
# gt3 = gct.GCT()
# gt3.build_from_DataFrame(cT)

#####################################
### plot DOS heatmap - make hist  ###
#####################################

cUnstack = cliqDos.unstack()
plt.hist(cUnstack,30)
plt.ylabel('freq',fontweight='bold')
plt.xlabel('median ' + matrixType,fontweight='bold')
plt.title('distribution of DOS-clique connections')
outF = os.path.join(wkdir, 'dos_clique_connection_distribution.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

### heatmap 
ccol.set_color_map()
plt.imshow(cliqDos.values,
    interpolation='nearest',
    aspect='auto')
# tickRange = range(0,tmpClust.shape[0])
# xtcks = [x for x in tmpClust.index]
# plt.xticks(tickRange, xtcks,rotation=90)
# plt.yticks(np.arange(len(xtcks)),xtcks)
plt.colorbar()
outF = os.path.join(wkdir, 'dos_clique_heatmap.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

############################################
### Examine cliq connection with members ###
############################################

n_members = 300
cpSer = pd.Series([item for sublist in cliqFrm['sig'] for item in sublist])
cpSer = cpSer[cpSer.isin(cliqFull.index)]
iRand = np.random.choice(len(cpLst),n_members,replace=False)
cpRand = np.random.choice(cpSer,n_members,replace=False)
rCliq = cliqFull.reindex(index=cpRand)


