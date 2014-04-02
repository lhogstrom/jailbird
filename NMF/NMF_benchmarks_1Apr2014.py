'''
-examine NMF results across cell lines
-build benchmarks

Larson Hogstrom, 12/2013
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.util.progress as update
from matplotlib import cm
from matplotlib.patches import Polygon

wkdir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation2/NMF_benchmark_development'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)
sourceDir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation2'

# directory of NMF result prefix and matrix dimentions
# dimDict = {'LINCS_core_c9_LM':'n4716x978',
# 'LINCS_core_c9_bing':'n4713x10638',
# dimDict = { 'PC3_c20_LM':'n585x978',
# 'PC3_c20_INF':'n585x10638',
# 'MCF7_c20_INF':'n652x10638',
# 'MCF7_c20_LM':'n652x978',
# 'MCF7_c9_INF':'n652x10638',
# 'MCF7_c9_LM':'n652x978',
# 'PC3_c9_INF':'n585x10638',
# 'PC3_c9_LM':'n585x978'}    

dimDict = {'A375_c9_lm_epsilon':'n473x978',
'A549_c9_lm_epsilon':'n612x978', # 
'HA1E_c9_lm_epsilon':'n578x978',
'HCC515_c9_lm_epsilon':'n543x978',
'HEPG2_c9_lm_epsilon':'n357x978',
'HT29_c9_lm_epsilon':'n433x978',
'MCF7_c9_lm_epsilon':'n652x978',
'PC3_c9_lm_epsilon':'n585x978',
'VCAP_c9_lm_epsilon':'n574x978'}

# dimDict = {'A375_c20_lm_epsilon':'n473x978',
# 'A549_c20_lm_epsilon':'n612x978', # 
# 'HA1E_c20_lm_epsilon':'n578x978',
# 'HCC515_c20_lm_epsilon':'n543x978',
# 'HEPG2_c20_lm_epsilon':'n357x978',
# 'HT29_c20_lm_epsilon':'n433x978',
# 'MCF7_c20_lm_epsilon':'n652x978',
# 'PC3_c20_lm_epsilon':'n585x978',
# 'VCAP_c20_lm_epsilon':'n574x978'}


sigDict = {} # significance counts
for prefix in dimDict:
    print prefix
    dim = dimDict[prefix]
    # local dir 
    graphDir = wkdir + '/' + prefix
    if not os.path.exists(graphDir):
        os.mkdir(graphDir)
    ### Load W and H matrix ###
    nComponents = 9
    Hfile = sourceDir + '/' + prefix + '/clique_compound_classes_' + dim + '.H.k' + str(nComponents) + '.gct'
    WFile = sourceDir + '/' + prefix + '/clique_compound_classes_' + dim + '.W.k' + str(nComponents) + '.gct'
    aFile = sourceDir + '/' + prefix + '/clique_compound_classes.v2.txt'
    Hmtrx = pd.read_csv(Hfile,sep='\t',skiprows=[0,1],index_col=0) #,header=True
    Hmtrx = Hmtrx.drop('Description',1)
    Hmtrx = Hmtrx.T
    anntFrm = pd.read_csv(aFile,sep='\t',index_col=0,header=None)
    headers= ['cell','cc','ss','is_gold','group_id','group','pert_id','group_name','tp','sig2']
    anntFrm.columns = headers
    anntFrm.index.name = 'sig1'
    # drop extra rows
    anntFrm = anntFrm[anntFrm.index.isin(Hmtrx.index)] # leave out annotations not in matrix
    ### read in mutual information matrices
    mFile = sourceDir + '/' + prefix + '/clique_compound_classes.MI.input_space.gct'
    mi = pd.read_csv(mFile,sep='\t',skiprows=[0,1],index_col=0) #,header=True
    mi = mi.drop('Description',1)
    cFile = sourceDir + '/' + prefix + '/clique_compound_classes.MI.k' + str(nComponents) + '.gct'
    cmi = pd.read_csv(cFile,sep='\t',skiprows=[0,1],index_col=0) #,header=True
    cmi = cmi.drop('Description',1)
    ### load in clique annotations and matrix
    cFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140221/cliques.gmt'
    cliqueGMT = gmt.read(cFile)
    cliqFrm = pd.DataFrame(cliqueGMT)
    #########################################
    ### graph individual group components ###
    #########################################
    maxVal = Hmtrx.max(axis=1).max()
    for r in cliqFrm.iterrows():
        grp = r[1]['id']
        brds = r[1]['sig']
        anntMtch = anntFrm[anntFrm.pert_id.isin(brds)]
        grpH = Hmtrx.reindex(anntMtch.index)
        meanVec = grpH.describe().ix['mean']
        ### take top three components - order acording to their strenght
        iTop3 = meanVec.order(ascending=False).index[:3]
        sortedTop = grpH.ix[:,iTop3].sort()
        topSum = sortedTop.sum(axis=1).order(ascending=False)
        grpH = grpH.ix[topSum.index,:] # sort acording to corr with mean
        Hfloat = np.float64(grpH.values)
        fig = plt.figure(figsize=(20, 10), dpi=50)
        plt.imshow(Hfloat,
            interpolation='nearest',
            cmap=cm.gray_r,
            vmax=maxVal)
        ytcks = list(grpH.index)
        xtcks = list(grpH.columns)
        plt.xticks(np.arange(len(xtcks)), xtcks,rotation=75)
        plt.yticks(np.arange(len(ytcks)),ytcks)
        plt.colorbar()
        plt.title(grp + ' - NMF component weights')
        grpMod = grpMod = ''.join(e for e in grp if e.isalnum())
        outF = os.path.join(graphDir,grpMod+'.png')
        plt.savefig(outF, bbox_inches='tight')
        plt.close()
    ##############################
    ### top component analysis ###
    ##############################
    # take the mean of the top 3 components for each group member
    topMeanDict = {}
    for r in cliqFrm.iterrows():
        grp = r[1]['id']
        brds = r[1]['sig']
        anntMtch = anntFrm[anntFrm.pert_id.isin(brds)]
        grpH = Hmtrx.reindex(anntMtch.index)
        meanVec = grpH.describe().ix['mean']
        #get top 
        nTop = 3 # number of top largest components to sort by
        iTop3 = meanVec.order(ascending=False).index[:nTop]
        sortedTop = grpH.ix[:,iTop3].sort()
        topSum = sortedTop.sum(axis=1).order(ascending=False)
        topMeanDict[grp] = topSum.mean()
    topMeanSer = pd.Series(topMeanDict)
    ##############################
    ### build null distribution ##
    ##############################
    # shuffle signatures frtopom random drugs - keep same group size
    nPerm = 4000
    zFrm = np.zeros([cliqFrm.shape[0],nPerm])
    nullMean = pd.DataFrame(zFrm,index=cliqFrm['desc'])
    prog = update.DeterminateProgressBar('cliq group')
    for irr,r in enumerate(cliqFrm.iterrows()):
        grp = r[1]['id']
        prog.update(grp,irr,len(cliqFrm.desc))
        brds = r[1]['sig']
        anntMtch = anntFrm[anntFrm.pert_id.isin(brds)]
        for ir in range(nPerm):
            nGrp = anntMtch.shape[0]
            iRand = np.random.choice(Hmtrx.index.values,nGrp,replace=False)
            grpH = Hmtrx.reindex(iRand)
            meanVec = grpH.mean()
            #get mean of top components
            nTop = 3 # number of top largest components to sort by
            iTop3 = meanVec.order(ascending=False).index[:nTop]
            sortedTop = grpH.ix[:,iTop3].sort()
            topSum = sortedTop.sum(axis=1).order(ascending=False)
            nullMean.ix[grp,ir] = topSum.mean()
    ##############################
    ### calculate significance  ##
    ##############################
    #compare each observed score to the null distribution
    pvalDict = {}
    for r in cliqFrm.iterrows():
        grp = r[1]['id']
        brds = r[1]['sig']
        obs = topMeanSer[grp]
        rndVec = nullMean.ix[grp,:]
        pvalDict[grp] = sum(rndVec > obs) / float(nPerm)
    pvalSer = pd.Series(pvalDict)
    pvalSer.name = 'top3_group_component_means'
    pvalSer.index.name = 'drug_group'
    pvalSer.sort()
    outF = graphDir + '/top_3_components_mean_group_pvalue.txt'
    pvalSer.to_csv(outF,sep='\t',header=True)
    # graph p-values
    fig = plt.figure(1, figsize=(14, 10))
    plt.plot(pvalSer,'.')
    outF = graphDir + '/top_3_components_mean_group_pvalues.png'
    xtcks = list(pvalSer.index.values)
    plt.xticks(np.arange(len(xtcks)), xtcks,rotation=90)
    plt.ylabel('p-value')
    plt.xlabel('pharmalogical class')
    plt.title(prefix + ' intra-class NMF component consistency')
    plt.savefig(outF, bbox_inches='tight')
    plt.close()
    ###
    sigDict[grp] = sum(pvalSer < .05) # number of groups w/ p-value under .05
    ##############################
    ### mutual information heatmaps ##
    ##############################
    # input space
    n = cmi.shape[0]
    mtrx = mi.ix[:n,:n]
    plt.imshow(mtrx,
        interpolation='nearest',
        cmap=cm.RdBu,
        vmin=-1,
        vmax=1)
    plt.xlabel('CMAP signatures')
    plt.ylabel('CMAP signatures')
    plt.colorbar()
    plt.title('mutual information - input space')
    outF = os.path.join(graphDir,'MI_matrix_input_space.png')
    plt.savefig(outF, bbox_inches='tight')
    plt.close()
    # nmf components
    n = cmi.shape[0]
    mtrx = cmi.ix[:n,:n]
    plt.imshow(mtrx,
        interpolation='nearest',
        cmap=cm.RdBu,
        vmin=-1,
        vmax=1)
    plt.xlabel('CMAP signatures')
    plt.ylabel('CMAP signatures')
    plt.colorbar()
    plt.title('mutual information - NMF components')
    outF = os.path.join(graphDir,'MI_matrix_NMF_components.png')
    plt.savefig(outF, bbox_inches='tight')
    plt.close()
    ##############################
    ### pairwise sigature comparisons ##
    ##############################
    # pairwise connections from input space
    inMI, outMI = MI_pairwise_comp(mi,cliqFrm,anntFrm)
    # pairwise connections from NMF component space
    inCMI, outCMI = MI_pairwise_comp(cmi,cliqFrm,anntFrm)
    ##############################
    ### simple boxplot ##
    ##############################
    ### NMF components - simple boxplot of intra-group connections
    fig = plt.figure(figsize=(8, 10), dpi=50)
    plt.boxplot(inCMI,vert=0)
    plt.xlim((-1,1))
    tickList = [x for x in inCMI.index]
    plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
    plt.tick_params(labelsize=8)
    plt.xlabel('mutual information',fontweight='bold')
    plt.title('intra-group connection - NMF components',fontweight='bold')
    outF = os.path.join(graphDir,'pairwise_comparison_boxplot_NMF_components.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()
    ### input space - simple boxplot of intra-group connections
    fig = plt.figure(figsize=(8, 10), dpi=50)
    plt.boxplot(inMI,vert=0)
    plt.xlim((-1,1))
    tickList = [x for x in inMI.index]
    plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
    plt.tick_params(labelsize=8)
    plt.xlabel('mutual information',fontweight='bold')
    plt.title('intra-group connection - input space',fontweight='bold')
    outF = os.path.join(graphDir,'pairwise_comparison_boxplot_input_space.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()    
    ##############################
    ### boxplot with null ##
    ##############################
    # alternate in-out arrays in a list 
    combSer = pd.Series()
    for ix,x in enumerate(inSorted):
        grp = inSorted.index[ix]
        combSer.set_value(grp+'_internal', inSorted[ix].values)
        combSer.set_value(grp+'_external', outSorted[ix].values)
    ### complex boxplot
    fig, ax1 = plt.subplots(figsize=(15,15))
    bp = plt.boxplot(combSer, notch=0, sym='+', vert=0, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='.')
    boxColors = ['darkkhaki','royalblue']
    numBoxes = len(combSer)
    medians = range(numBoxes)
    for i in range(numBoxes):
      box = bp['boxes'][i]
      boxX = []
      boxY = []
      for j in range(5):
          boxX.append(box.get_xdata()[j])
          boxY.append(box.get_ydata()[j])
      boxCoords = zip(boxX,boxY)
      # Alternate between Dark Khaki and Royal Blue
      k = i % 2
      boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
      ax1.add_patch(boxPolygon)
      # Now draw the median lines back over what we just filled in
      med = bp['medians'][i]
      medianX = []
      medianY = []
      for j in range(2):
          medianX.append(med.get_xdata()[j])
          medianY.append(med.get_ydata()[j])
          plt.plot(medianX, medianY, 'k')
          medians[i] = medianY[0]
    tickList = [x for x in combSer.index]
    plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
    plt.tick_params(labelsize=8)
    # plt.ylabel('compound class',fontweight='bold')
    plt.xlabel('mutual information',fontweight='bold')
    plt.title('intra-group connection - input space',fontweight='bold')
    outF = os.path.join(graphDir,'pairwise_comparison_boxplot_NMF_components_inter_intra.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()
    ##############################
    ### internal-external connection metric ##
    ##############################
    # intra-connection mean - inter-connection mean
    diffDict = {}
    for grp in inSorted.index:
        iMean = np.mean(inSorted[grp])
        oMean = np.mean(outSorted[grp])
        diff = iMean-oMean
        diffDict[grp] = diff   
    diffSer = pd.Series(diffDict)
    diffSer = diffSer[inSorted.index] # reorder acording to median connection
    # make plot
    plt.plot(diffSer, '.')
    tickList = [x for x in diffSer.index]
    plt.xticks(np.arange(1,len(tickList)+1),tickList,rotation=90)
    plt.tick_params(labelsize=8)
    plt.ylabel('difference in mutual information',fontweight='bold')
    # plt.title('intra-group connection - input space',fontweight='bold')
    outF = os.path.join(graphDir,'mutual_information_difference_metric.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()
    # save to text file
    outF = os.path.join(graphDir,'mutual_information_difference_metric.txt')
    diffSer.name = 'intra_inter_group_MI_difference'
    diffSer.index.name = 'group'
    diffSer.to_csv(outF, sep='\t',header=True)



def MI_pairwise_comp(MI_matrix,cliqFrm,anntFrm):
    '''
    -take loop through compound groups 
    -sort median mutual information within group connections
    -sort median mutual information of group members to non-members

    Parameters
    ----------
    MI_matrix : pandas dataFrame
        pairwise mutual information of signatures
    cliqFrm : pandas dataFrame
        drug class asignments
    anntFrm : pandas dataFrame
        signature annotation

    returns
    ----------
    inSorted : pandas Series
        -intra group pairwise mutual information scores
        -sorted by intr-group MI median
    outSorted : pandas Series
        -group members to non-members pairwise mutual information
        -sorted by intr-group MI median

    ''' 
    inList = {} # store intra-group mutual information
    outList = {} # store inter-group mutual information
    for r in cliqFrm.iterrows():
        grp = r[1]['id']
        brds = r[1]['sig']
        anntMtch = anntFrm[anntFrm.pert_id.isin(brds)]
        anntNonMtch = anntFrm[~anntFrm.pert_id.isin(brds)]
        # within group comparisons
        grpIn = MI_matrix.reindex(index=anntMtch.index, columns=anntMtch.index)
        ilRand = np.triu_indices(len(grpIn),k=0)
        upIn = grpIn.values.copy()
        upIn[ilRand] = np.nan
        vGI = upIn[~np.isnan(upIn)]
        inList[grp] = vGI
        # outside group comparisons 
        grpOut = MI_matrix.reindex(index=anntNonMtch.index, columns=anntMtch.index)
        vGO = grpOut.unstack()
        outList[grp] = vGO.values
    # inMedian = [np.median(x) for x in inList]
    inSer = pd.Series(inList)
    inMedian = inSer.apply(np.median)
    inMedian.sort()
    inMedian = inMedian[~np.isnan(inMedian)]
    # inSorted = inSer[inMedian.index].values
    inSorted = inSer[inMedian.index]
    outSer = pd.Series(outList)
    outSorted = outSer[inMedian.index]
    return (inSorted, outSorted)
