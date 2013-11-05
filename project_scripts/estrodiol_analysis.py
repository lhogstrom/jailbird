'''
analyze estrodiol compounds 

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

# get directory
dir1 = '/xchip/cogs/projects/pharm_class' 
wkdir = dir1 + '/estrodiol_analysis_Oct11'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

pDescDict = {'BRD-A60070924':'17-alpha-estradiol',
    'BRD-K04218075':'clomiphene',
    'BRD-K95309561':'dienestrol',
    'BRD-K45330754':'diethylstilbestrol',
    'BRD-K04046242':'equilin',
    'BRD-A74907996':'equol',
    'BRD-K18910433':'estradiol-17-beta',
    'BRD-K17016787':'estriol',
    'BRD-K81839095':'estrone',
    'BRD-A83237092':'fulvestrant',
    'BRD-K43797669':'genistein',
    'BRD-K63828191':'raloxifene',
    'BRD-K93754473':'tamoxifen',
    'BRD-K67174588':'toremifene'}
estrodiolBrds = pDescDict.keys()

#
for brd in estrodiolBrds:
    out = wkdir + '/' pDescDict[brd]
    if not os.path.exists(out):
        os.mkdir(out)
    # examine drug rosiglitazone, on CPC006 plates, gold signatures
    # query = {'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True} 
    query = {'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True, 'pert_dose':{'$gt':3}}
    self = PertExplorer(pert_id = brd,
                        metric = 'wtcs',
                        query = query,
                        out = out)
    # cluster the data, store it on the PertExplorer instance
    self.hclust = HClust(self.data.copy(), self.annots.copy(), out = self.out,
                         symmfun = np.max, cutoff = 0.5)
    self.hclust.cluster()
    self.hclust.get_cluster_order()
    self.hclust.draw_dendrogram(orientation = 'right')
    self.hclust.save_dendrogram()
    clust_order = self.hclust.cluster_order.order().index
    # add cluster assignment to annotations
    self.annots['hclust_assignment'] = self.hclust.cluster_assignment
    # draw heatmap, ordered by the clustering; save
    self.draw_heatmap(order = clust_order,
                      vmin = 0, vmax = 1,
                      cmap = cm.OrRd,
                      col_label = 'cell_id',
                      row_label = None,
                      annots_top = [('hclust_assignment', False)],
                      annots_right = ['pert_dose',
                                      ('distil_cc_q75', {'vmin' : 0, 'vmax' : 1}),
                                      ('distil_ss', {'vmin' : 0, 'vmax' : 10})],
                      annots_bottom = [('cell_lineage', True)],
                      showfig = False,
                      title = 'Example clustering for sig_id ' + pDescDict[brd])
    self.heatmap.save(format = 'png')
    # draw baseline data; save
    self.showBaselines(order = self.annots.cell_id[clust_order].values,
                       clust_assignment = self.hclust.cluster_assignment.order().values)
    self.saveBaselines(format = 'png')
    # write html file and annotations
    self.write_html(extra_figs = [(self.hclust.dendro_fname, 'dendrogram')])
    self.write_annotations(encoding = 'latin-1')

### run sig_introspect on each compound seperatly
max_processes=7
processes = set()
sigIdDict = {}
# make S-C plots, color by cell line
for brd in estrodiolBrds:
    # out = wkdir + '/' + pDescDict[brd]
    # if not os.path.exists(out):
    #     os.mkdir(out)
    CM = mu.CMapMongo()
    #all cell lines
    # qRes = CM.find({'pert_id':brd,'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True},
    #         {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
    #         toDataFrame=True)
    #MCF7 only query
    qRes = CM.find({'pert_id':brd,'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True,'cell_id':'MCF7'},
            {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
            toDataFrame=True)
    if type(qRes) == int:
        continue
    sigIDs = list(qRes['sig_id'].values)
    sigIdDict[brd] = sigIDs
    # make SC plots
    out = wkdir + '/SC_plots_MCF7'
    if not os.path.exists(out):
        os.mkdir(out)
    sco = sc.SC()
    sco.add_sc_from_mongo_sig_ids(sigIDs)
    scOut = out + '/' + pDescDict[brd] + '_sc.png'
    sco.plot(out=scOut)
    ### run sig_intorspect
    outIntrospect = out + '/brew_introspect_results_MCF7'
    # outIntrospect = out + '/brew_introspect_results_all_cell_lines'
    if not os.path.exists(outIntrospect):
        os.mkdir(outIntrospect)
    qSer = qRes['sig_id']
    iname = pDescDict[brd]
    outF = outIntrospect + '/' + iname + '_sig_ids.grp'
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
    # processes.add(subprocess.Popen(cmd,shell=True))
    # if len(processes) >= max_processes:
    #     os.wait()
    #     processes.difference_update(
    #         p for p in processes if p.poll() is not None)

### run sig_introspect on all compounds together
cellList = ['HEPG2', 'A375', 'PC3', 'MCF7', 'HA1E', 'HT29', 'HCC515', 'A549']
nonMCF7cells = ['HEPG2', 'A375', 'PC3', 'HA1E', 'HT29', 'HCC515', 'A549']
outIntrospect = wkdir + '/brew_introspect_nonMCF7'
if not os.path.exists(outIntrospect):
    os.mkdir(outIntrospect)
CM = mu.CMapMongo()
qRes = CM.find({'pert_id':{'$in': estrodiolBrds},'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True,'cell_id':{'$in': nonMCF7cells}},
        {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
        toDataFrame=True)
qSer = qRes['sig_id']
outF = outIntrospect + '/nonMCF7_sig_ids.grp'
qSer.to_csv(outF,index=False,header=False)
#run sig_introspect
cmd = ' '.join(['rum -q local',
     '-d sulfur_io=100',
     '-o ' + outIntrospect,
     '-x sig_introspect_tool ',
     '--sig_id ' + outF,
     '--metric wtcs',
     '--out ' + outIntrospect])
os.system(cmd)

# sim_set - pert_mfc
cellList = ['HEPG2', 'A375', 'PC3', 'MCF7', 'HA1E', 'HT29', 'HCC515', 'A549']
for cell in cellList:
    CM = mu.CMapMongo()
    qRes = CM.find({'pert_id':{'$in': estrodiolBrds},'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True,'cell_id':cell},
            {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
            toDataFrame=True)
    # look at pert_id-pert_mfc_id pairs
    pairs = qRes['pert_id'] + ':' + qRes['pert_mfc_id']
    pairSet = set(pairs.values)
    ### run sig_introspect on all of the sig_ids at once
    outIntrospect = wkdir + '/brew_introspect_all_cps_' + cell
    # outIntrospect = out + '/brew_introspect_results_all_cell_lines'
    if not os.path.exists(outIntrospect):
        os.mkdir(outIntrospect)
    qSer = qRes['sig_id']
    outF = outIntrospect + '/all_cps_' + cell + '_sig_ids.grp'
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


# get instances from  sig_info - distil_id
CM = mu.CMapMongo()
sigList = qRes['sig_id']
sigInfo = mu.sig_info(sigList)
distilDict = {}
for sig in sigInfo:
    distilDict[sig['sig_id']] = sig['distil_id']
# # loop through each compound and write grp of roast ids
# for brd in estrodiolBrds:
#     iname = pDescDict[brd]
#     # cp dir
#     out = wkdir + '/' + iname
#     if not os.path.exists(out):
#         os.mkdir(out)
#     # introspect dir
#     outIntrospect = out + '/roast_introspect_MCF7'
#     if not os.path.exists(outIntrospect):
#         os.mkdir(outIntrospect)        
#     if brd in sigIdDict:
#         sigs = sigIdDict[brd]
#         distilIDList = []
#         for sig in sigs:
#             distilIDList.extend(distilDict[sig])
#         dIDSer = pd.Series(distilIDList)
#         outF = outIntrospect + '/' + iname + '_distil_ids.grp'
#         dIDSer.to_csv(outF,index=False,header=False)
#         #run sig_introspect - with roast data
#         cmd = ' '.join(['rum -q hour',
#              '-d sulfur_io=100',
#              '-o ' + outIntrospect,
#              '-x sig_introspect_tool',
#              '--build_id custom',
#              '--score ' + cmap.roast_path,
#              '--cid ' + outF,
#              '--metric wtcs',
#              '--out ' + outIntrospect])
#         os.system(cmd)

### put distil_ids into gmt form --> each row is a different (sig_id? or pert_csf_id?)

### make list of dictionaries for making gmt file
CM = mu.CMapMongo()
sigIdDict = {}
for brd in estrodiolBrds:
    #all cell lines
    qRes = CM.find({'pert_id':brd,'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True},
            {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
            toDataFrame=True)
    #MCF7 only query
    # qRes = CM.find({'pert_id':brd,'sig_id' : {'$regex' : 'PCLB00'}, 'is_gold' : True,'cell_id':'MCF7'},
    #         {'sig_id':True,'pert_iname':True,'pert_id':True,'pert_mfc_id':True},
    #         toDataFrame=True)
    if type(qRes) == int:
        continue
    sigIDs = list(qRes['sig_id'].values)
    sigIdDict[brd] = sigIDs
gmtList = []
for brd in sigIdDict:
    print brd
    sigs = sigIdDict[brd]
    distilIDList = []
    for sig in sigs:
        distilIDList.extend(distilDict[sig])
    # dIDSer = pd.Series(distilIDList)
    gmtDict = {}
    gmtDict['id'] = brd
    gmtDict['desc'] = pDescDict[brd]
    gmtDict['sig'] = distilIDList
    gmtList.append(gmtDict)
gmtOut = wkdir + '/roast_ids_all_cell_lines.gmt'
gmt.write(gmtList,gmtOut)


# look up ESR1 and ESR2 baseline expression in each cell line
MC = mu.MongoContainer()
ci = MC.gene_info.find({'pr_gene_symbol':'ESR2'},{},toDataFrame=True)
# exprDict = MC.gene_info.find({'pr_gene_symbol':'ESR1','pr_id':'205225_at'},{'is_expressed':True})
exprDict = MC.gene_info.find({'pr_gene_symbol':'ESR2','pr_id':'210780_at'},{'is_expressed':True})
exprSer = pd.Series(exprDict[0])
exprSer.ix[cellList]

# affyDict = MC.gene_info.find({'pr_gene_symbol':'ESR1','pr_id':'205225_at'},{'basex_affx':True})
affyDict = MC.gene_info.find({'pr_gene_symbol':'ESR2','pr_id':'210780_at'},{'basex_affx':True})
affySer = pd.Series(affyDict[0])
affySer.ix[cellList]

# rSeqDict = MC.gene_info.find({'pr_gene_symbol':'ESR1','pr_id':'205225_at'},{'basex_rnaseq':True})
rSeqDict = MC.gene_info.find({'pr_gene_symbol':'ESR2','pr_id':'210780_at'},{'basex_rnaseq':True})
rSeqSer = pd.Series(rSeqDict[0])
rSeqSer.ix[cellList]

### load in results from sig_introspect - all cps, all cell liens
ispecDir = '/xchip/cogs/projects/pharm_class/estrodiol_analysis_Oct10/brew_introspect_all_cps_all_cell_lines/oct10/my_analysis.sig_introspect_tool.2013101010482180'
rnkptFile = '/self_rankpt_n266x266.gctx'
inF = ispecDir + rnkptFile
rnkptGCT = gct.GCT()
rnkptGCT.read(inF)
rnkpt = rnkptGCT.frame
# organize cells and sig_ids
sigIds = rnkpt.index.values
cells = [x.split('_')[1] for x in sigIds]
brdColum = [x.split(':')[1] for x in sigIds]
inames = [pDescDict[x[:13]] for x in brdColum]
cellSer = pd.Series(sigIds,index=cells)
sigSer = pd.Series(sigIds,index=sigIds)
sigFrm = pd.DataFrame(sigSer,columns=['sig_id'])
sigFrm['brd'] = brdColum
sigFrm['cell'] = cells 
sigFrm['iname'] = inames
# within cell connections
MCF7vals = []
notMCF7vals = []
tickList = []
rnkptSumList = []
for cell in set(cells):
    cellSigs = cellSer.ix[cell]
    cellMtrx = rnkpt.ix[cellSigs,cellSigs]
    upMtrx = av_mtrx_upper(cellMtrx.values)
    flatMtrx = upMtrx.flatten()
    rnkptArray = flatMtrx[~np.isnan(flatMtrx)]
    rnkptSumList.append(rnkptArray)
    if cell == 'MCF7':
        MCF7vals.extend(rnkptArray)
    else:
        notMCF7vals.extend(rnkptArray)
    tickList.append(cell)
### within compound, within cell line
grped = sigFrm.groupby(['brd','cell'])
GrpSigDict = grped.groups
keys = GrpSigDict.keys()
for cell in set(cells):
    cellKeys = [x for x in keys if x[1] == cell]
    tickList = []
    rnkptSumList = []
    for key in cellKeys:
        iname = pDescDict[key[0][:13]]
        sigs = GrpSigDict[key]
        cellMtrx = rnkpt.ix[sigs,sigs]
        upMtrx = av_mtrx_upper(cellMtrx.values)
        flatMtrx = upMtrx.flatten()
        rnkptArray = flatMtrx[~np.isnan(flatMtrx)]
        if len(rnkptArray) > 0:
            rnkptSumList.append(rnkptArray)
            tickList.append(iname)
    plt.boxplot(rnkptSumList,vert=0)
    plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
    plt.tick_params(labelsize=8)
    plt.ylabel('compound',fontweight='bold')
    plt.xlabel('rank point',fontweight='bold')
    plt.title(cell + ' self connection of estriodiol and related compounds',fontweight='bold')
    outF = os.path.join(wkdir, cell + '_estriodiol_anologs_boxplot.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()

### plot raw values
bMin = 150
bMax = 200
binwidth = 1
tmpBins = range(bMin,bMax+binwidth,binwidth)
Histbins = np.divide(tmpBins,float(2))
plt.hist(notMCF7vals,bins=Histbins,normed=True,alpha=.3,color='b')
plt.hist(MCF7vals,bins=Histbins,normed=True,alpha=.3,color='r')
plt.xlim([70,100])
plt.show()
### plot square
sqMCF7 = np.power(MCF7vals,7)
sqNmcf7 = np.power(notMCF7vals,7)
plt.hist(sqNmcf7,30,normed=True,alpha=.3,color='b')
plt.hist(sqMCF7,30,normed=True,alpha=.3,color='r')
# plt.xlim([0,np.power(1000,4)])
plt.show()





#sumscore boxplot
plt.boxplot(rnkptSumList,vert=0)
plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
plt.tick_params(labelsize=8)
plt.ylabel('cell line',fontweight='bold')
plt.xlabel('rank point',fontweight='bold')
plt.title('connection of estriodiol and related compounds',fontweight='bold')
outF = os.path.join(wkdir, 'estriodiol_anologs_boxplot.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)


def av_mtrx_upper(mtrx):
    '''
    take the average of a pairwise matrix - return upper right of matrix
    
    '''        
    nm = len(mtrx)
    avMtrx = np.zeros((nm,nm))
    for i1 in range(nm):
        for i2 in range(nm):
            val1 = mtrx[i1,i2]
            val2 = mtrx[i2,i1]
            if (val1 == -666) or (val2 == -666):
                avMtrx[i1,i2] = np.nan
            else:
                avMtrx[i1,i2] = stats.nanmean([val1,val2])            
    # avMtrxUp = np.triu(avMtrx,k=1)
    iUp = np.tril_indices(nm)
    avMtrx[iUp] = np.nan
    return avMtrx



### make a different boxplot for each cell line (cp self connections)