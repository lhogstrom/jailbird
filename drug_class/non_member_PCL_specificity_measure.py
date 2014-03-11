'''
-identify which compounds do not already belong to PCL groups
-which non-member compounds show strong, specific connections to PCL groups

Larson Hogstrom, 2/2014
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
import cmap.plot.colors as ccol
import scipy.cluster
import  HTML

wkdir = '/xchip/cogs/projects/pharm_class/lhwork/non_member_analysis_5Mar2014'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

########################
### retrieve cp info ###
########################

### use lass matrix
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_lass_n39560x7147.gctx'
# matrixType = 'rnkpt_indp_lass'
### use mrp4 mtrx
# mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/indep_mrp4_n39560x7147.gctx'
# matrixType = 'mrp4'
### use lass matched matrix
mtrxSummly = '/xchip/cogs/projects/connectivity/summly/matrices/matched_lass_n7147x7147.gctx'
matrixType = 'rnkpt_matched_lass'
# nullMtrx = '/xchip/cogs/projects/connectivity/null/dmso/lass_n1000x7147.gctx'
nullMtrx = '/xchip/cogs/projects/connectivity/null/random/lass_n1000x7147.gctx'

##########################################
### get brds for the whole matrix ###
##########################################

gt = gct.GCT()
gt.read(mtrxSummly)
summFrm = gt.frame
pInames = gt.get_column_meta('pert_iname')
pIDs = gt.get_column_meta('pert_id')
pType = gt.get_column_meta('pert_type')
# anntFrm = pd.DataFrame({'pert_id':pIDs,'pert_type':pType,'pert_iname':pInames},index=pIDs)
pInameType = []
for i,x in enumerate(pInames):
    pInameType.append(pInames[i]+ '.' +pType[i])
anntFrm = pd.DataFrame({'pert_id':pIDs,'pert_type':pType,'pert_iname':pInames},index=pInameType)
sigSer = pd.Series(index=summFrm.index, data=summFrm.columns)
outGRP = wkdir + '/summly_matched_ids.grp'
sigSer.to_csv(outGRP,index=False)

####################
### load cliques ###
####################

# groupGMT = '/xchip/cogs/projects/pharm_class/rnwork/cliques/cpd_groups_n147.gmt'
groupGMT = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140213/cliques.gmt'
# groupGMT = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140221/cliques.gmt'

cliqueGMT = gmt.read(groupGMT)
cliqFrm = pd.DataFrame(cliqueGMT)
cliqFrm['group_size'] = cliqFrm.sig.apply(len)
cliqFrm.index = cliqFrm['desc']
 cliqFrm['Name'].str.replace("/","-")

# which compounds are clique members vs. non-members
cliqMemberLong = [item for sublist in cliqFrm.sig.values for item in sublist]
cliqMemb = list(set(cliqMemberLong))
isMemb = anntFrm.pert_id.isin(cliqMemb)
isCp = anntFrm.pert_type == 'trt_cp'
nonMemb = anntFrm[isCp & ~isMemb].index.values
nonMid = anntFrm[isCp & ~isMemb].pert_id.values

#########################################
### load sig_cliquescore_tool results ###
#########################################

# cliques against all compounds
# cFileFull = '/xchip/cogs/projects/DOS/PCL_comparison_Feb162014/feb14/my_analysis.sig_cliquescore_tool.2014021411101091/clique_median_n145x7147.gctx'
cFileFull = '/xchip/cogs/projects/DOS/PCL_comparison_Feb142014/feb14/my_analysis.sig_cliquescore_tool.2014021414592391/clique_median_n80x7147.gctx'
gt2 = gct.GCT()
gt2.read(cFileFull)
cliqFull = gt2.frame
cliqFull = cliqFull.reindex(sigSer.values)
cliqFull.index = sigSer.index
# transpose matrix and write
# cT = cliqFull.T 
# gt3 = gct.GCT()
# gt3.build_from_DataFrame(cT)

# cliques against DMSO
# cFile = '/xchip/cogs/projects/DOS/PCL_comparison_Feb142014/feb14/my_analysis.sig_cliquescore_tool.2014021415010091/clique_median_n80x1000.gctx'
# gt4 = gct.GCT()
# gt4.read(cFile)
# cliqDMSO = gt4.frame

# # cliques against random signatures
# cFile = '/xchip/cogs/projects/DOS/PCL_comparison_Feb142014/feb14/my_analysis.sig_cliquescore_tool.2014021415003091/clique_median_n80x1000.gctx'
# gt5 = gct.GCT()
# gt5.read(cFile)
# cliqRnd = gt5.frame

#########################################
### sig_cliquescore_tool heatmap ###
#########################################

# nonMGroup = [item for sublist in nonMembDict.values() for item in sublist]
# nmCliq = cliqFull.reindex(matchFrm.index)
nmCliq = cliqFull.reindex(nonMid)
# make heatmap
plt.close()
ccol.set_color_map()
fig = plt.figure(1, figsize=(20, 25))
plt.imshow(nmCliq.T.values,
    interpolation='nearest',
    aspect='auto')
xtickRange = range(0,nmCliq.shape[0])
xtcks = [x for x in nmCliq.index]
ytickRange = range(0,nmCliq.shape[1])
ytcks = [x for x in nmCliq.columns]
# plt.xticks(xtickRange, xtcks,rotation=90)
plt.yticks(ytickRange, ytcks)
plt.xlabel('compounds')
plt.title('median connection of non-member cps to cliques')
plt.colorbar()
ctFile = os.path.join(wkdir, 'non_member_cp_clique_heatmap.png')
plt.savefig(ctFile, bbox_inches='tight',dpi=200)
plt.close()

#########################################
### make specificity index for sig_cliquescore_tool  ###
#########################################

# cliquescore results for non-clique members
cliqNM = cliqFull.reindex(nonMid)
cThresh = 85
cPass = cliqNM > cThresh
cpSum = cPass.sum(axis=1) # number of clique connections above threshold
plt.hist(cpSum)
plt.ylabel('# cliq connections',fontweight='bold')
plt.xlabel('compound',fontweight='bold')
plt.title('clique connections above lass score ' + str(cThresh))
outF = os.path.join(wkdir, 'compound_cliq_connection_sum.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

# specificity ratio: highest score: next highest 5 scores
cliqPass = cliqNM[cpSum > 0]
# cliqPass = cliqNM[cpSum == 1]
def specif_ratio(x):
    x = x.order(ascending=False)
    max1 = x[0]
    # max2 = x[1:5] 
    max2 = x[1] 
    topTwoRatio = max1/max2
    pHigher = (max1-max2)/max2
    #ratio:
    # highest value
    # median of next 5 highest?
    return pHigher
ratioSer = cliqPass.apply(specif_ratio,axis=1)

plt.hist(ratioSer,30)
plt.ylabel('ratio of highest cliq connections',fontweight='bold')
plt.xlabel('compound',fontweight='bold')
plt.title('specificity of clique connections above lass score ' + str(cThresh))
outF = os.path.join(wkdir, 'compound_cliq_connection_specificity.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()
highRatio = ratioSer[ratioSer > .5]
highRatio.sort(ascending=False)
highSpecificity = cliqPass.ix[highRatio.index,:]

#########################################
### heatmap of high specificity clique hits ###
#########################################

#rows
mtrx = highSpecificity
Y = scipy.cluster.hierarchy.linkage(mtrx, method='centroid')
Z = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
cOrder = Z['leaves']
iPCL = mtrx.index[cOrder]
clustered = mtrx.reindex(index=iPCL)
#columns 
Y_ = scipy.cluster.hierarchy.linkage(mtrx.T, method='centroid')
Z_ = scipy.cluster.hierarchy.dendrogram(Y,orientation='right')
cOrder_ = Z['leaves']
iPCL_columns = mtrx.columns[cOrder]
clustered = mtrx.reindex(index=iPCL,columns=iPCL_columns)
#heatmap
plt.close()
ccol.set_color_map()
fig = plt.figure(1, figsize=(20, 25))
plt.imshow(mtrx.T.values,
    interpolation='nearest',
    aspect='auto')
xtickRange = range(0,mtrx.shape[0])
xtcks = [x for x in mtrx.index]
ytickRange = range(0,mtrx.shape[1])
ytcks = [x for x in mtrx.columns]
plt.xticks(xtickRange, xtcks,rotation=90)
plt.yticks(ytickRange, ytcks)
plt.xlabel('high specificity compounds')
plt.title('median connection of non-member cps to cliques')
plt.colorbar()
ctFile = os.path.join(wkdir, 'non_member_cp_clique_high_specificity.png')
plt.savefig(ctFile, bbox_inches='tight',dpi=200)
plt.close()

#####################################
### make Clique heatmaps for hits ###
#####################################
#top hit dictionary
matchDict = {}
for cpTup in highSpecificity.iterrows():
    cp = cpTup[0]
    cpSer = cpTup[1]
    matchDict[cp] = cpSer.idxmax()
# heatmaps
summFrm.columns = summFrm.index
imageDict = {}
for cp,cliq in matchDict.iteritems():
    # add compound to larger list of brds
    grpCps = cliqFrm.ix[cliq,:].sig
    grpCps = grpCps[:]
    grpCps.append(cp)
    grpFrm = summFrm.reindex(index=grpCps,columns=grpCps)
    # grpAnnt = anntFrm.reindex(grpCps)
    grpAnnt = anntFrm[anntFrm.pert_id.isin(grpCps)]
    grpAnnt.index = grpAnnt.pert_id
    grpAnnt = grpAnnt.reindex(grpCps)
    grpFrm.index = grpAnnt.pert_iname
    cpIname = grpAnnt.pert_iname[-1]
    # grpFrm.index = grpAnnt.pert_iname
    # score between outside member to member
    pairVec = grpFrm.ix[:,cp][:-1]
    medSumm = pairVec.median()
    matchDict[cp] = cliq
    ### make heatmap
    plt.close()
    ccol.set_color_map()
    fig = plt.figure(1, figsize=(5, 5))
    plt.imshow(grpFrm.values,
        interpolation='nearest',
        aspect='auto',
        vmin=-100,
        vmax=100)
    xtickRange = range(0,grpFrm.shape[0])
    xtcks = [x for x in grpFrm.index]
    plt.xticks(xtickRange, xtcks,rotation=90)
    plt.yticks(xtickRange, xtcks)
    plt.xlabel('compounds')
    plt.title(cliq + ' compounds plus - ' + cp + ' (' + cpIname+ ')')
    plt.colorbar()
    FilePrefix = cliq + '_' + cp + '_potential_member_heatmap.png'
    FilePrefix = FilePrefix.replace("/","-")
    outF = os.path.join(wkdir, FilePrefix)
    imageDict[cp] = FilePrefix
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close()

##############################
### make heatmap sub-pages ###
##############################

imageSer = pd.Series(imageDict)
imageSer.name = 'heatmap_image'
htmlDict = {}
for x in imageSer.iteritems():
    brd = x[0]
    img = x[1]
    prefix = img[:-3] + 'html'
    pageF = os.path.join(wkdir,prefix)
    htmlDict[brd] = prefix
    with open(pageF,'w') as f:
            lineWrite = '<BR>' + '<tr><td><img src=' + img + '></td></tr> <BR>'
            f.write(lineWrite + '\n')
prefixSer = pd.Series(htmlDict)
htmlSer = '<a href=./' + prefixSer + '>heatmap</a> <BR>'
htmlSer.name = 'heatmap_image'
# make table 
matchSer = pd.Series(matchDict)
matchSer.name = 'PCL_hit_group'
matchFrm = pd.DataFrame(matchSer)
anntMtch = anntFrm[anntFrm.pert_id.isin(matchFrm.index)]
anntMtch.index = anntMtch.pert_id
matchFrm = pd.concat([matchFrm,anntMtch],axis=1)
matchFrm = pd.concat([matchFrm,pd.DataFrame(htmlSer)],axis=1)
matchFrm = matchFrm.sort('PCL_hit_group')
hitFile = wkdir + '/summly_dendrogram_non_member_hits.txt'
matchFrm.to_csv(hitFile,sep='\t',index=False)
# write html page
pageF = os.path.join(wkdir,'index.html')
with open(pageF,'w') as f:
    ## write hit table
    table_data = []
    with open(hitFile,'rt') as Fread:
        for line in Fread:
            table_data.append(line.split('\t'))
    htmlcodeFDR = HTML.table(table_data)
    f.write(htmlcodeFDR + '\n')
    # write heatmap src
    lineWrite = '<BR>' + '<tr><td><img src="non_member_cp_clique_high_specificity.png" height="1500" width="1000"></td></tr> <BR>'
    f.write(lineWrite + '\n')
