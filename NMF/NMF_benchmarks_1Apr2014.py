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
import matplotlib.gridspec as gridspec

########################
### define functions ###
########################

def MI_pairwise_comp(MI_matrix,cliqFrm,anntFrm,graphDir):
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

def boxplot_with_null(inSorted,outSorted,graphDir,space_name='input_space'):
    '''
    -make boxplot of intra-group connections compare these to 
    how the group signatures connect with non-group signatures

    Parameters
    ----------
    inSorted : pandas Series
        -intra group pairwise mutual information scores
        -sorted by intr-group MI median
    outSorted : pandas Series
        -group members to non-members pairwise mutual information
        -sorted by intr-group MI median
    ''' 
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
    plt.title('intra-group connection - ' + space_name,fontweight='bold')
    outF = os.path.join(graphDir,'pairwise_comparison_boxplot_NMF_components_inter_intra_' + space_name + '.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()

def intra_group_boxplot(inSorted,graphDir,space_name='input_space'):
    '''
    -make boxplot of intra-group connections 

    Parameters
    ----------
    inSorted : pandas Series
        -intra group pairwise mutual information scores
        -sorted by intr-group MI median
    outSorted : pandas Series
        -group members to non-members pairwise mutual information
        -sorted by intr-group MI median
    ''' 
    fig = plt.figure(figsize=(8, 10), dpi=50)
    plt.boxplot(inSorted,vert=0)
    plt.xlim((-1,1))
    tickList = [x for x in inSorted.index]
    plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
    plt.tick_params(labelsize=8)
    plt.xlabel('mutual information',fontweight='bold')
    plt.title('intra-group connection - ' + space_name,fontweight='bold')
    outF = os.path.join(graphDir,'pairwise_comparison_boxplot_ ' + space_name + '.png')
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
    plt.title('intra-group connection - ' + space_name,fontweight='bold')
    outF = os.path.join(graphDir,'pairwise_comparison_ ' + space_name + '.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()    

def MI_pairwise_metric(inSorted,outSorted,graphDir,space_name='input_space'):
    '''
    -take loop through compound groups 
    -calculate difference in mean mutual information from within group connections
    compared to mean mutual information of group members to non-members

    Parameters
    ----------
    inSorted : pandas Series
        -intra group pairwise mutual information scores
        -sorted by intr-group MI median
    outSorted : pandas Series
        -group members to non-members pairwise mutual information
        -sorted by intr-group MI median

    returns
    ----------
    diffSer : pandas Series
        difference in mean mutual information from within group connections
        compared to mean mutual information of group members to non-members

    ''' 
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
    plt.ylim((-.3,1))
    tickList = [x for x in diffSer.index]
    plt.xticks(np.arange(1,len(tickList)+1),tickList,rotation=90)
    plt.tick_params(labelsize=8)
    plt.ylabel('difference in mutual information',fontweight='bold')
    plt.title('intra-group connection difference - ' + space_name ,fontweight='bold')
    outF = os.path.join(graphDir,'mutual_information_difference_metric_' + space_name + '.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()
    # save to text file
    outF = os.path.join(graphDir,'mutual_information_difference_metric_' + space_name + '.txt')
    diffSer.name = 'intra_inter_group_MI_difference'
    diffSer.index.name = 'group'
    diffSer.to_csv(outF, sep='\t',header=True)
    return diffSer

def mi_heatmap(mtrx,graphDir,space_name='input_space'):
    plt.imshow(mtrx,
        interpolation='nearest',
        cmap=cm.RdBu,
        vmin=-1,
        vmax=1)
    plt.xlabel('CMAP signatures')
    plt.ylabel('CMAP signatures')
    plt.colorbar()
    plt.title('mutual information - ' + space_name)
    outF = os.path.join(graphDir,'MI_matrix_' + space_name + '.png')
    plt.savefig(outF, bbox_inches='tight')
    plt.close()

def mi_heatmap_group(mtrx,graphDir,grpSort,space_name='input_space'):
    fig = plt.figure(1, figsize=(10, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[40, 1]) 
    ax0 = plt.subplot(gs[0])
    ax0.imshow(mtrx,
        interpolation='nearest',
        cmap=cm.RdBu_r,
        vmin=-1,
        vmax=1,
        aspect='auto')
    plt.xlabel('CMAP signatures')
    plt.ylabel('CMAP signatures')
    plt.title('mutual information - ' + space_name)
    smFrm = pd.DataFrame(grpSort.group_id.astype(np.float))
    ax1 = plt.subplot(gs[1])
    ax1.imshow(smFrm,
        interpolation='nearest',
        cmap=cm.jet,
        aspect='auto')
    ax1.axes.get_yaxis().set_visible(False)
    plt.xticks([0], ['drug_group'],rotation=90)
    outF = os.path.join(graphDir,'MI_matrix_' + space_name + '.png')
    plt.savefig(outF, bbox_inches='tight')
    plt.close()

def group_component_maps(Hmtrx,cliqFrm,graphDir):
    '''
    -for all signatures of a given class, make a heatmap of the NMF
    H-matrix weights
    -order acording to the mean of the top 3 most heavily weighted components

    Parameters
    ----------
    Hmtrx : pandas DataFrame
        -matrix of NMF weightings for each signatures (n_signatures x n_components)
    cliqFrm : pandas DataFrame
        -signature annotationsw
    '''
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

def combine_group_top_components(Hmtrx,cliqFrm,nTop=3,metric='mean'):
    '''
    -take the mean of the top 3 components for members of each group

    Parameters
    ----------
    Hmtrx : pandas DataFrame
        -matrix of NMF weightings for each signatures (n_signatures x n_components)
    cliqFrm : pandas DataFrame
        -signature annotationsw
    nTop : int
        -number of top largest components to sort by

    Returns
    ----------
    topMeanFrm : pandas DataFrame
        for each group:
            group_signature_counts  
            top_mean_metric
    '''
    topMeanDict = {}
    sigCountDict = {}
    for r in cliqFrm.iterrows():
        grp = r[1]['id']
        brds = r[1]['sig']
        anntMtch = anntFrm[anntFrm.pert_id.isin(brds)]
        grpH = Hmtrx.reindex(anntMtch.index)
        if grpH.shape[0] < 2: # skip groups with two or less signatures
            continue
        else:
            sigCountDict[grp] = grpH.shape[0]
            meanVec = grpH.describe().ix[metric]         
            iTop3 = meanVec.order(ascending=False).index[:nTop]
            sortedTop = grpH.ix[:,iTop3].sort()
            topSum = sortedTop.sum(axis=1).order(ascending=False)
            topMeanDict[grp] = topSum.mean()
    topMeanFrm = pd.DataFrame({'top_mean_metric':topMeanDict,'group_signature_counts':sigCountDict})
    return topMeanFrm

def build_combine_null(Hmtrx,cliqFrm,topMeanFrm,nTop=3,nPerm=4000):
    '''
    -
    shuffle signatures from random drugs - keep same group size

    Parameters
    ----------
    Hmtrx : pandas DataFrame
        -matrix of NMF weightings for each signatures (n_signatures x n_components)
    cliqFrm : pandas DataFrame
        -signature annotationsw
    nTop : int
        -number of top largest components to sort by  
    topMeanFrm : pandas DataFrame
        for each group:
            group_signature_counts  
            top_mean_metric         

    Returns
    ----------
    nullMean : pandas DataFrame
        -matrix of null distributions for each group size
    '''
    # count the number of signatures in the input data that belong to each group
    cliqSize = topMeanFrm.group_signature_counts
    groupSizeSet = set(cliqSize) # unique group sizes across input groups
    zFrm = np.zeros([len(groupSizeSet),nPerm])
    nullMean = pd.DataFrame(zFrm,index=np.sort(list(groupSizeSet))) # one row for each group size
    nullMean.index.name = 'group_size'
    prog = update.DeterminateProgressBar('group size')
    for ix,nGrp in enumerate(nullMean.index):
        prog.update(nGrp,ix,len(nullMean.index))
        for ir in range(nPerm):
            iRand = np.random.choice(Hmtrx.index.values,nGrp,replace=False)
            grpH = Hmtrx.reindex(iRand)
            meanVec = grpH.mean()
            #get mean of top components
            iTop3 = meanVec.order(ascending=False).index[:nTop]
            sortedTop = grpH.ix[:,iTop3].sort()
            topSum = sortedTop.sum(axis=1).order(ascending=False)
            nullMean.ix[nGrp,ir] = topSum.mean()
    return nullMean

def top_components_combined_significance(topMeanFrm,nullMean,nPerm=4000):
    '''
    -compare each observed score to the null distribution
    -make graph of pvalues found for each group

    Parameters
    ----------
    topMeanFrm : pandas DataFrame
        for each group:
            group_signature_counts  
            top_mean_metric
    nullMean : pandas DataFrame
        -matrix of null distributions for each group size
    nPerm : int
        number of permutations used in 'build_combine_null'

    Returns
    ----------
    pvalSer : pandas Series
        -vector of p-values
    '''
    pvalDict = {}
    for grp in topMeanFrm.index:
        obs = topMeanFrm.top_mean_metric[grp]
        nGrp = topMeanFrm.group_signature_counts[grp]
        rndVec = nullMean.ix[nGrp,:]
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
    return pvalSer

def group_rolling_density(cliqFrm,anntFrm,cluster_mtrx,window=10,group_min=4,space_name='input_space'):
    '''
    -for each group, explore a clustered input matrix
    -evaluate the density of group membership along the clustered axis
    -identify group 'hotspots' within a given window

    Parameters
    ----------
    cliqFrm : pandas dataFrame
        drug class asignments
    anntFrm : pandas dataFrame
        signature annotation

    cluster_mtrx : pandas DataFrame
        -matrix of clustered signatures
    window : int
        density window in which to count group members
    group_min : int
        minimum number of group members within a window in which to call 
        significant

    Returns
    ----------
    cluster_density_max : pandas Series
        -maximum 'hotspot' size for group members within a given window size
    '''
    densityDict = {}
    # prog = update.DeterminateProgressBar('cliq graph')
    for ir,r in enumerate(cliqFrm.iterrows()):
        grp = r[1]['id']
        brds = r[1]['sig']
        # prog.update(grp,ir,len(cliqFrm.index))
        anntMtch = anntFrm[anntFrm.pert_id.isin(brds)]
        iSigs = anntMtch.index
        nSigs = anntMtch.shape
        if nSigs > 2: # skip groups with two or less signatures
            boolAr = cluster_mtrx.index.isin(iSigs)
            boolSer = pd.Series(index=cluster_mtrx.index, data=boolAr)
            rollSum = pd.stats.moments.rolling_sum(boolSer,window)
            rollSum.name = 'rolling_group_count'
            nrollSum = rollSum[~rollSum.isnull()]
            rollFrm = pd.DataFrame(rollSum)
            rollFrm['location'] = np.arange(len(rollSum))
            rollFrm[grp] = boolSer
            densityDict[grp] = max(nrollSum)
            if max(nrollSum) > group_min:
                # print grp
                peakSer = nrollSum[nrollSum > group_min]
                rollFrm.ix[nrollSum[nrollSum > group_min].index,:]
                rMin = rollFrm.index.get_loc(peakSer.index[0])
                rMax = rollFrm.index.get_loc(peakSer.index[-1])
                rollRange = np.arange(rMin,rMax)
                if len(rollRange) <= 100: # skip if cluster is too long
                    localFrm = rollFrm.ix[(rMin-10):(rMax+10),:]
                    out = graphDir + '/' + grp + '_local_dendrogram_table_' + space_name + '.txt'
                    localFrm.to_csv(out,sep='\t')
                # create graph
                out = graphDir + '/' + grp + '_cluster_rolling_sum_' + space_name + '.png'
                plt.plot(rollSum)
                plt.ylim((0,window))
                plt.xlabel('cluster axis')
                plt.ylabel('rolling sum')   
                plt.title(grp + ' - target density from clustering - window = ' +str(window)) 
                plt.savefig(out, bbox_inches='tight')
                plt.close()
    cluster_density_max = pd.Series(densityDict)
    cluster_density_max.sort(ascending=False)
    cluster_density_max.name = 'group_cluster_density_max'
    cluster_density_max.index.name = 'group'
    outF = graphDir + '/cluster_rolling_density_max_' + space_name + '.txt'
    cluster_density_max.to_csv(outF,sep='\t')
    ### graph max values
    fig = plt.figure(1, figsize=(10, 10))
    plt.plot(cluster_density_max,'.')
    outF = graphDir + '/cluster_rolling_density_max_' + space_name + '.png'
    xtcks = list(cluster_density_max.index.values)
    plt.xticks(np.arange(len(xtcks)), xtcks,rotation=90)
    plt.ylabel('max group density')
    plt.xlabel('pharmalogical class')
    plt.title('group density in clustered connection matrix - ' + space_name)
    plt.savefig(outF, bbox_inches='tight')
    plt.close()
    return cluster_density_max

##################
### set inputs ###
##################

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

# dimDict = {'A375_c9_lm_epsilon':'n473x978',
# 'A549_c9_lm_epsilon':'n612x978', # 
# 'HA1E_c9_lm_epsilon':'n578x978',
# 'HCC515_c9_lm_epsilon':'n543x978',
# 'HEPG2_c9_lm_epsilon':'n357x978',
# 'HT29_c9_lm_epsilon':'n433x978',
# 'MCF7_c9_lm_epsilon':'n652x978',
# 'PC3_c9_lm_epsilon':'n585x978',
# 'VCAP_c9_lm_epsilon':'n574x978'}

dimDict = {'A375_c20_lm_epsilon':'n473x978',
'A549_c20_lm_epsilon':'n612x978', # 
'HA1E_c20_lm_epsilon':'n578x978',
'HCC515_c20_lm_epsilon':'n543x978',
'HEPG2_c20_lm_epsilon':'n357x978',
'HT29_c20_lm_epsilon':'n433x978',
'MCF7_c20_lm_epsilon':'n652x978',
'PC3_c20_lm_epsilon':'n585x978',
'VCAP_c20_lm_epsilon':'n574x978'}

sigDict = {} # significance counts
for prefix in dimDict:
    print prefix
    dim = dimDict[prefix]
    # local dir 
    graphDir = wkdir + '/' + prefix
    if not os.path.exists(graphDir):
        os.mkdir(graphDir)
    ### Load W and H matrix ###
    nComponents = 20
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
    # group_component_maps(Hmtrx,cliqFrm,graphDir)
    # # ##############################
    # # ### top component analysis ###
    # # ##############################
    # # take the mean of the top 3 components for each group member
    # topMeanFrm = combine_group_top_components(Hmtrx,cliqFrm,metric='mean')
    # # ##############################
    # # ### build null distribution ##
    # # ##############################
    # # repeate metric - but shuffle signatures from groups of equal size
    # nullMean = build_combine_null(Hmtrx,cliqFrm,topMeanFrm,nTop=3,nPerm=4000)
    # ##############################
    # ### calculate significance  ##
    # ##############################
    # pvalSer = top_components_combined_significance(topMeanFrm,nullMean,nPerm=4000)
    # sigDict[grp] = sum(pvalSer < .05) # number of groups w/ p-value under .05
    # ##############################
    # ### mutual information heatmaps ##
    # ##############################
    # input space
    # mi_heatmap(mi,graphDir,space_name='input_space')
    # mi_heatmap(cmi,graphDir,space_name=str(nComponents) +'_NMF_components')
    # ### reindex matrix by group
    # grpSort = anntFrm.sort('group')
    # # by input
    # miGrp = mi.ix[grpSort.index, grpSort.index]
    # # mi_heatmap(miGrp,graphDir,space_name='group_ordered_input_space')
    # mi_heatmap_group(miGrp,graphDir,grpSort,space_name='group_ordered_input_space')
    # # components
    # cmiGrp = cmi.ix[grpSort.index, grpSort.index]
    # # mi_heatmap(cmiGrp,graphDir,space_name=str(nComponents) +'_NMF_components_by_group')
    # mi_heatmap_group(cmiGrp,graphDir,grpSort,space_name=str(nComponents) +'_NMF_components_by_group')
    # ##############################
    # ### pairwise sigature comparisons ##
    # ##############################
    # # pairwise connections from input space
    inMI, outMI = MI_pairwise_comp(mi,cliqFrm,anntFrm,graphDir)
    # pairwise connections from NMF component space
    inCMI, outCMI = MI_pairwise_comp(cmi,cliqFrm,anntFrm,graphDir)
    # ##############################
    # ### simple boxplot ##
    # ##############################
    # ### NMF components - simple boxplot of intra-group connections

    ##############################
    ### boxplot with null ##
    ##############################
    # alternate in-out arrays in a list 
    boxplot_with_null(inMI,outMI,graphDir,space_name='input_space')
    boxplot_with_null(inCMI,outCMI,graphDir,space_name=str(nComponents) +'_NMF_components_by_group')
    ##############################
    ### internal-external connection metric ##
    ##############################
    # intra-connection mean - inter-connection mean
    # diffSer = MI_pairwise_metric(inCMI,outCMI,graphDir,space_name='NMF_components')
    # diffSer = MI_pairwise_metric(inMI,outMI,graphDir,space_name='input_space')
    ##############################
    ### rolling sum of group density in clustered matrix ##
    ##############################
    # cluster_density_max = group_rolling_density(cliqFrm,anntFrm,mi,window=10,group_min=7,space_name='input_space')
    # cluster_density_max = group_rolling_density(cliqFrm,anntFrm,cmi,window=10,group_min=7,space_name=str(nComponents) +'_NMF_components_by_group')



# query benchmark 
# - what percent of top 20 connections is a group-member
# - what percent of the group's signatures are in the top 50?

# clustering benchmark
# - clustering density
