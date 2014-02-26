'''
-identify which compounds do not already belong to PCL groups
-look at the summly dendrogram, which compounds apear close to the cp cluster
-how specfic is the non-member relationship to to clique members?

Larson Hogstrom, 2/2014
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
import scipy.cluster
import cmap.util.progress as progress

wkdir = '/xchip/cogs/projects/pharm_class/lhwork/non_member_analysis_25Feb2014'
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
cliqueGMT = gmt.read(groupGMT)
cliqFrm = pd.DataFrame(cliqueGMT)
cliqFrm['group_size'] = cliqFrm.sig.apply(len)
cliqFrm.index = cliqFrm['desc']

# which compounds are clique members vs. non-members
cliqMemberLong = [item for sublist in cliqFrm.sig.values for item in sublist]
cliqMemb = list(set(cliqMemberLong))
isMemb = anntFrm.pert_id.isin(cliqMemb)
isCp = anntFrm.pert_type == 'trt_cp'
nonMemb = anntFrm[isCp & ~isMemb].index.values

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

###############################
### load summly tree #### 
###############################

### load XML file 
sTreeXML = '/xchip/cogs/projects/connectivity/clustering/matched_lass_dendro.xml'
tstPanLst = open(sTreeXML).read().splitlines()
labelList = []
for str1 in tstPanLst:
    if 'label="' in str1:
        s1 = str1.split('label="')[1]
        s2 = s1.split('">')[0]
        labelList.append(s2)
labelList.pop(0)
clustFrm = anntFrm.reindex(labelList)
clustFrm['order'] = range(clustFrm.shape[0])
outF = wkdir + '/summly_dendrogram_table.txt'
clustFrm.to_csv(outF,sep='\t')


##########################################
### rolling clique sum along dendrogram ###
##########################################

### rolling sum window
graphDir = wkdir + '/window_min'
if not os.path.exists(graphDir):
    os.mkdir(graphDir)
window=10
group_min = 4 # minimum clique members within a window
clustered_groups = {}
nonMembDict = {}
prog = progress.DeterminateProgressBar('cliq graph')
for icliq,cliq in enumerate(cliqFrm.desc):
    prog.update(cliq,icliq,len(cliqFrm.desc))
    cliqMod = cliqFrm.ix[icliq,'desc']
    brds = cliqFrm.ix[icliq,'sig']
    boolSer = clustFrm.pert_id.isin(brds)
    isNonMemb = clustFrm.index.isin(nonMemb)
    rollSum = pd.stats.moments.rolling_sum(boolSer,window)
    rollSum.name = 'rolling_group_count'
    nrollSum = rollSum[~rollSum.isnull()]
    rollFrm = pd.DataFrame(rollSum)
    rollFrm['location'] = np.arange(len(rollSum))
    rollFrm[cliqMod] = boolSer
    rollFrm['is_non_clique_member'] = isNonMemb
    rollFrm['pert_id'] = clustFrm.pert_id
    if max(nrollSum) > group_min:
        peakSer = nrollSum[nrollSum > group_min]
        clustered_groups[cliq] = peakSer
        rollFrm.ix[nrollSum[nrollSum > group_min].index,:]
        rMin = rollFrm.index.get_loc(peakSer.index[0])
        rMax = rollFrm.index.get_loc(peakSer.index[-1])
        rollRange = np.arange(rMin,rMax)
        if len(rollRange) <= 100: # skip if cluster is too long
            localFrm = rollFrm.ix[(rMin-10):(rMax+10),:]
            out = graphDir + '/' + cliqMod + '_local_dendrogram_table.txt'
            localFrm.to_csv(out,sep='\t')
            # any local dos compounds
            if localFrm.is_non_clique_member.any():
                nonMembDict[cliq] = list(localFrm[localFrm.is_non_clique_member].pert_id)
                nonMembDict[cliq] = list(localFrm[localFrm.is_non_clique_member].pert_id)
        # create graph
        out = graphDir + '/' + cliqMod + '_rolling_sum.png'
        plt.plot(rollSum)
        plt.ylim((0,window))
        plt.xlabel('cluster axis')
        plt.ylabel('rolling sum')   
        plt.title(cliq + ' - target density from clustering - window = ' +str(window)) 
        plt.savefig(out, bbox_inches='tight')
        plt.close()

# list of DOS compounds that are close by clique region of the dendrogram
dosDendro = pd.Series(nonMembDict)
outF = wkdir + '/non_clique_member_cps_proximal_to_local_cliques.txt'
dosDendro.to_csv(outF,sep='\t')

#########################################
### dendrogram analysis ###
#########################################


###
matchDict = {}
imageDict = {}
summFrm.columns = summFrm.index
prog = progress.DeterminateProgressBar('cliq graph')
for ix,x in enumerate(nonMembDict.iteritems()):
    cliq = x[0]
    cps = x[1]
    prog.update(cliq,ix,len(nonMembDict))
    for cp in cps:
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
        if medSumm > 90:
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
            outF = os.path.join(graphDir, FilePrefix)
            imageDict[cp] = FilePrefix
            plt.savefig(outF, bbox_inches='tight',dpi=200)
            plt.close()
### make heatmap sub-pages
imageSer = pd.Series(imageDict)
imageSer.name = 'heatmap_image'
htmlDict = {}
for x in imageSer.iteritems():
    brd = x[0]
    img = x[1]
    prefix = img[:-3] + 'html'
    pageF = os.path.join(graphDir,prefix)
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
pageF = os.path.join(graphDir,'index.html')
with open(pageF,'w') as f:
    ## write hit table
    table_data = []
    with open(hitFile,'rt') as Fread:
        for line in Fread:
            table_data.append(line.split('\t'))
    htmlcodeFDR = HTML.table(table_data)
    f.write(htmlcodeFDR + '\n')
    # lineWrite = '<BR>' + '<tr><td><img src=./non_member_cp_clique_heatmap.png></td></tr> <BR>'
    # f.write(lineWrite + '\n')
    # for imageFile in fileList:
    #     # csF = os.path.join(outdir,cell1 + '_drug_target_CS_dist.png')
    #     lineWrite = '<BR>' + '<tr><td><img src=' + imageFile + '></td></tr>'
    #     f.write(lineWrite + '\n')


#########################################
### sig_cliquescore_tool heatmap ###
#########################################

# nonMGroup = [item for sublist in nonMembDict.values() for item in sublist]
nmCliq = cliqFull.reindex(matchFrm.index)
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
plt.xticks(xtickRange, xtcks,rotation=90)
plt.yticks(ytickRange, ytcks)
plt.xlabel('compounds')
plt.title('median connection of non-member cps to cliques')
plt.colorbar()
ctFile = os.path.join(graphDir, 'non_member_cp_clique_heatmap.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

