'''
analyze TRIB1 compounds 

including these tasks:
1) count sig_id and experimental paramaters for four compounds
2) make SC plots for each compound
3) look up modification of TRIB1 LM gene from compounds
4) pairwise summly comparison 
5) summly indpenedent mode

October 2013
'''
import numpy as np
from cmap.analytics.pert_explorer import PertExplorer
from cmap.analytics.cluster import HClust
import cmap.analytics.sc as sc
import cmap
import os
from os import path
from matplotlib import cm
import cmap.util.mongo_utils as mu
import subprocess
import cmap.tools.sig_dose_tool as sdt
import cmap.io.gct as gct
import pandas as pd
import cmap.io.gmt as gmt
from scipy import stats
import matplotlib.pyplot as plt
import cmap.plot.colors as colors
import cmap

# get directory
dir1 = '/xchip/cogs/projects/TRIB1' 
wkdir = dir1 + '/TRIB1_analysis_Oct21'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

#define compounds of interest
trib1Cps = ['BRD-K75627148', 
    'BRD-K35860134', 
    'BRD-K67774729', 
    'BRD-K16956545', 
    'BRD-K16410418']

### 1 ) get signature info
CM = mu.CMapMongo()
qRes = CM.find({'pert_id':{'$in':trib1Cps}, 'is_gold' : True},
        {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True,'cell_id':True,'pert_time':True,'is_gold':True},
        toDataFrame=True)
qRes.index = qRes['sig_id']
outF = wkdir + '/TRIB1_signature_details.txt'
qRes.to_csv(outF,sep='\t',index=False)

grped = qRes.groupby('pert_id')
for grp in grped.groups:
    nSig = len(grped.groups[grp])
    # print grp + ' ' + str(nSig)
    print str(nSig)

### descriptive data on signatures
#number of cell lines per compound
grped = qRes.group_by('pert_iname')


### 2) make S-C plots, color by cell line
sigIdDict = {}
for brd in trib1Cps:
    # out = wkdir + '/' + pDescDict[brd]
    # if not os.patmh.exists(out):
    #     os.mkdir(out)
    CM = mu.CMapMongo()
    #all cell lines
    qRes = CM.find({'pert_id':brd, 'is_gold' : True,},
            {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True,'cell_id':True,'pert_time':True},
            toDataFrame=True)
    if type(qRes) == int:
        continue
    sigIDs = list(qRes['sig_id'].values)
    sigIdDict[brd] = sigIDs
    # make SC plots
    out = wkdir + '/SC_plots'
    if not os.path.exists(out):
        os.mkdir(out)
    sco = sc.SC()
    sco.add_sc_from_mongo_sig_ids(sigIDs)
    scOut = out + '/' + brd + '_sc.png'
    sco.plot(out=scOut,title=brd)

### 3) look up expression modulation of TRIB1
# what is the affy id for TRIB1
mc = mu.MongoContainer()
geneInf = mc.gene_info.find({'pr_gene_symbol':'TRIB1'},
        {'pr_gene_symbol':True,'pr_id':True},
        toDataFrame=True)
trib1_affy_ID = geneInf['pr_id'][0]
sigs = qRes['sig_id'].values
# look up Z scores
afPath = cmap.score_path
gt = gct.GCT()
gt.read(src=afPath,cid=sigs,rid='lm_epsilon')
zFrm = gt.frame
trib1Z = zFrm.ix[trib1_affy_ID,:]
trib1Z.name = 'TRIB1_z_score'
#add new field to table
qRes = qRes.merge(pd.DataFrame(trib1Z),left_index=True,right_index=True)
# look up ranks
gtRank = gct.GCT()
gtRank.read(src=cmap.rank_lm_path,cid=sigs,rid='lm_epsilon')
rankFrm = gtRank.frame
trib1Rank = rankFrm.ix[trib1_affy_ID,:]
trib1Rank.name = 'TRIB1_rank'
qRes = qRes.merge(pd.DataFrame(trib1Rank),left_index=True,right_index=True)
# write table to file
outF = wkdir + '/TRIB1_signature_details.txt'
qRes.to_csv(outF,sep='\t',index=False)
grped = qRes.groupby(['pert_iname','cell_id'])
for grp in grped.groups:
    iGroup = grped.groups[grp]
    grpRes = qRes.ix[iGroup,:]
#look at distribution for each compound
for brd in trib1Cps:
#     # sigs = qRes['sig_id'].values
    sigRes = qRes[qRes['pert_id'] == brd]
    sigs = cellRes['sig_id'].values
    gt = gct.GCT()
    gt.read(src=afPath,cid=sigs,rid='lm_epsilon')
    zFrm = gt.frame
    trib1Z = zFrm.ix[trib1_affy_ID,:]
    ### make histogram of TRIB1 z score
    plt.hist(trib1Z,30)
    plt.xlim([-100,100])
    plt.ylabel('freq',fontweight='bold')
    plt.xlabel('z-score',fontweight='bold')
    plt.title('TRIB1 z scores - in response to ' + brd)
    outF = os.path.join(wkdir, brd + '_TRIB1_z_score_hist.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()
## trib rank

## other genes of interest - LDLR, PCSK9, MTTP

# t-test to compare TRIB1 regulation to null
tribVec = trib1Z.values
# nullVec = tribVec + 10
[tStat, pVal] = stats.ttest_ind(tribVec,nullVec)

## repeat in individual cell lines
# cellLines = set(qRes['cell_id'])
# for cell in cellLines:
#     # sigs = qRes['sig_id'].values
#     cellRes = qRes[qRes['cell_id'] == cell]

### 4) make pairwise summly matrix of trib1 cps
## load in summly matrix
rnkptGCT = gct.GCT()
# load merged matrix ranpoint matrix
rnkptGCT.read(cmap.summly_rankpt_path)
rnkpt = rnkptGCT.frame
rnkpt = rnkpt.replace(-666,np.nan)
cols = rnkpt.columns.values
colShort = [x.split('|')[0] for x in cols]
colShortSet = list(set(colShort))
rnkpt.columns = colShort #shorten column names
# index of first instance of each column
iColFirst = []
dupList = []
for ix,x in enumerate(colShort):
    if x not in dupList:
        iColFirst.append(ix)
        dupList.append(x)
rnkptNoDup = rnkpt.ix[:,iColFirst]
# full summly matrix
rnkptNoDup = rnkptNoDup.reindex(index=colShortSet,columns=colShortSet)
#summly matrix of trib1 compounds
testedMtrx = rnkptNoDup.reindex(index=trib1Cps,columns=trib1Cps)
#summly matrix - raw
outMtx = wkdir + '/TRIB1_summly_pairwise_mtrx_raw'
make_group_heatmap('TRIB1_compounds',testedMtrx.index.values, testedMtrx.values,outMtx)
#symetric summly matrix
avSummMtrx = av_mtrx(testedMtrx.values)
outMtx = wkdir + '/TRIB1_summly_pairwise_mtrx_symentric'
make_group_heatmap('TRIB1_compounds',testedMtrx.index.values, avSummMtrx,outMtx)

### test overlap among top connectors
nTop = 200
summRes = rnkptNoDup.reindex(columns=trib1Cps)
# topConnectors = pd.Series()
topConnectors = []
for brd in summRes.columns:
    print brd
    brdSer = summRes.ix[:,brd]
    brdOrder = brdSer.order(ascending=False)
    topList = brdOrder[:nTop]
    topConnectors.extend(topList.index.values)
    # topConnectors.append(topList)
connSer = pd.Series(topConnectors)
freqTop = connSer.value_counts()
topPerts = freqTop[freqTop>3].index.values
#get pert Inames for top Perts
cm = mu.CMapMongo()
topIname = cm.find({'pert_id':{'$in':list(topPerts)}},
                {'pert_id':True,'pert_iname':True},
                toDataFrame=True)
inameGrped = topIname.groupby('pert_id')
inameGrped.first()

### 5) summly independent mode


### 6) what are the other TRIB1 regulators in the DB?
# rank order perturbations in terms of TRIB1 regulation
mc = mu.MongoContainer()
geneInf = mc.gene_info.find({'pr_gene_symbol':'TRIB1'},
        {'pr_gene_symbol':True,'pr_id':True},
        toDataFrame=True)
trib1_affy_ID = geneInf['pr_id'][0]
sigs = qRes['sig_id'].values
# look up Z scores for TRIB1 
afPath = cmap.score_path
gt = gct.GCT(cmap.score_path)
rids = gt.get_gctx_rid()
iTrib1 = rids.index(trib1_affy_ID)
gt.read_gctx_matrix(row_inds=[iTrib1])
cids = gt.get_gctx_cid()
# make Frame/Series of TRIB1 regulation
trib1Ser = pd.Series(gt.matrix[0,:],index=cids)
trib1Ser.name = 'z_score'
#filter gold vs. not gold
CM = mu.CMapMongo()
goldQuery = CM.find({'is_gold' : True}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True},
        toDataFrame=True)
goldSigs = goldQuery['sig_id']
trib1Gold = trib1Ser.reindex(goldSigs.values)
# order and build Frame
trib1Sort = trib1Gold.order(ascending=False)
affSigs = trib1Sort.index.values
pert_ids = [x.split(':')[1] for x in affSigs]
pert_id_short = [x.split(':')[1][:13] for x in affSigs]
cell_ids = [x.split('_')[1] for x in affSigs]
trib1Frm = pd.DataFrame(trib1Sort)
trib1Frm['pert_id'] = pert_ids
trib1Frm['pert_id_short'] = pert_id_short
trib1Frm['cell_id'] = cell_ids
trib1Frm['rank'] = np.arange(len(trib1Frm)) + 1 
trib1Frm['percent_rank'] = trib1Frm['rank']/float(len(trib1Frm))
outF = wkdir + '/CMAP_top_regulating_signatures.txt'
topTrib1 = trib1Frm[trib1Frm['z_score'] >= 2]
topTrib1.to_csv(outF,sep='\t')

## where do the apriori compounds sit?
# matchFrm = trib1Frm[trib1Frm['pert_id_short'] == brd]
matchFrm = trib1Frm.reindex(qRes.index)
#plot overall distribution of TRIB1 connctions
plt.hist(trib1Sort,30,color='b',normed=True,alpha=.5,range=[-10,10],label='all signatures')
plt.hist(matchFrm['z_score'],30,color='r',normed=True,alpha=.5,range=[-10,10],label='qPCR hits')
plt.title('TRIB1 regulation')
plt.xlabel('TRIB1 z score')
plt.ylabel('normalized freq')
plt.legend()
outF = os.path.join(wkdir,'trib1_z_score_distribution.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()
# plot percent ranks of apriori compounds in TRIB1 reulation
plt.hist(matchFrm['percent_rank'],30,color='r',normed=False,alpha=.5,range=[0,1])
plt.title('5 compounds, 150 signatures of interest - TRIB1 percent rank \n relative to other is_gold signatures')
plt.xlabel('percent rank of TRIB1 regulation')
plt.ylabel('freq')
outF = os.path.join(wkdir,'trib1_regulation_percent_rank.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()


# group by cell line
largeGrped = trib1Frm.groupby('cell_id')
cell= 'HEPG2'
iCell = largeGrped.groups[cell]
cellFrm = trib1Frm.ix[iCell]
# recalculate the rankings for the specific cell line 
cellFrm['rank'] = np.arange(len(cellFrm)) + 1 
cellFrm['percent_rank'] = cellFrm['rank']/float(len(cellFrm))
cellOrder = cellFrm.sort('rank')
# cellOrder = cellOrder[cellOrder['z_score'] >= 2]
outF = wkdir + '/HEPG2_top_regulating_signatures.txt'
cellFrm.to_csv(outF,sep='\t')

#get results specific to a cell line
aprioriGrped = qRes.groupby('cell_id')
iCellSigs = aprioriGrped.groups[cell]
matchCellFrm = cellFrm.ix[iCellSigs]
#plot 
# binSet=np.logspace(-10, 0, 50)
plt.hist(cellFrm['z_score'],30,color='b',normed=True,alpha=.5,range=[-10,10],label='all signatures')
plt.hist(matchCellFrm['z_score'],30,color='r',normed=True,alpha=.5,range=[-10,10],label='qPCR hits')
plt.title('TRIB1 regulation in HEPG2')
plt.xlabel('TRIB1 z score')
plt.ylabel('normalized freq')
plt.legend()
outF = os.path.join(wkdir,'HEPG2_trib1_z_score_distribution.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()
# plot percent ranks of apriori compounds in TRIB1 reulation
plt.hist(matchCellFrm['percent_rank'],30,color='r',normed=False,alpha=.5,range=[0,1])
plt.title(cell + ' - 5 compounds, 33 signatures of interest \n TRIB1 percent rank relative to other is_gold signatures')
plt.xlabel('percent rank of TRIB1 regulation')
plt.ylabel('freq')
outF = os.path.join(wkdir,'HEPG2_trib1_regulation_percent_rank.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()


#put into the frame 

# get pert_inames and pert_ids
CM = mu.CMapMongo()
qRes = CM.find({'sig_id':{'$in':list(affSigs)}}, #, 'is_gold' : True
        {'sig_id':True,'pert_iname':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True},
        toDataFrame=True)

## double check the trib1Ser 
# gtCheck = gct.GCT()
# gtCheck.read(src=afPath,cid=sigs,rid=[trib1_affy_ID])
# gtCheck.frame.ix[:,'PCLB003_PC3_24H:BRD-K75627148-001-02-4:0.04']
# trib1Ser['PCLB003_PC3_24H:BRD-K75627148-001-02-4:0.04']



### 7) run sig_introspect results
### run sig_introspect on all compounds together
outIntrospect = wkdir + '/brew_introspect'
if not os.path.exists(outIntrospect):
    os.mkdir(outIntrospect)
CM = mu.CMapMongo()
qRes = CM.find({'pert_id':{'$in':trib1Cps}, 'is_gold' : True},
        {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True,'cell_id':True,'pert_time':True,'is_gold':True},
        toDataFrame=True)
qRes.index = qRes['sig_id']
qSer = qRes['sig_id']
outF = outIntrospect + '/TRIB1_cp_sig_ids.grp'
qSer.to_csv(outF,index=False,header=False)
#run sig_introspect
cmd = ' '.join(['rum -q hour',
     '-d sulfur_io=100',
     '-o ' + outIntrospect,
     '-x sig_introspect_tool ',
     '--sig_id ' + outF,
     '--metric wtcs',
     '--out ' + outIntrospect])
os.system(cmd)

def make_group_heatmap(group_name, grp,rankpt_mtrx, out):
    '''
    Write a matrix to file

    Parameters
    ----------
    group_name : str
        brds paired with inames
    grp : list
        compound brds that make up the group
    rankpt_mtrx : numpy.ndarray
        rankpt matrix
    sumRank_mtrx : numpy.ndarray
        matrix of percent summly rank scores
    out : str
        output path        
    ''' 
    fig = plt.figure(1, figsize=(20, 8))
    plt.suptitle(group_name + ' compound group',fontsize=14, fontweight='bold')
    plt.subplot(111)
    plt.title('mean_rankpt_4')
    colors.set_color_map()
    plt.imshow(rankpt_mtrx,
            interpolation='nearest',
            vmin=-100, 
            vmax=100)
    ytcks = [x for x in grp]
    plt.xticks(np.arange(len(grp)), ytcks,rotation=75)
    plt.yticks(np.arange(len(grp)),ytcks)
    plt.colorbar()
    fig.savefig(out, bbox_inches='tight')
    plt.close()

def av_mtrx(mtrx):
    '''
    take the average of a pairwise matrix
    
    '''        
    nm = len(mtrx)
    avMtrx = np.zeros((nm,nm))
    for i1 in range(nm):
        for i2 in range(nm):
            if i1 == i2:
                avMtrx[i1,i2] = np.nan
                continue
            val1 = mtrx[i1,i2]
            val2 = mtrx[i2,i1]
            if (val1 == -666) or (val2 == -666):
                avMtrx[i1,i2] = np.nan
            else:
                avMtrx[i1,i2] = stats.nanmean([val1,val2])
    return avMtrx
