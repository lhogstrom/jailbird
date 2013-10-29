'''
PCLA (pharmacological class analyzer)

Given a set of drugs, the tool examines the distribution of summly
connections scores within each drug class

Larson Hogstrom, 9/2013
'''

import numpy as np
import os
import cmap.util.mongo_utils as mu
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import matplotlib
from cmap.analytics.pert_explorer import PertExplorer
from cmap.analytics.cluster import HClust
import cmap.plot.colors as colors
import cmap.analytics.modzsig as modzsig
import cmap.io.gmt as gmt
import cmap.io.gct as gct
import cmap
from cmap.analytics.queryer import Queryer
from os import path
import glob 
import shutil
from scipy import stats

class PCLA(object):
    '''
    Class to serve as a container for pathway connections. 

    Parameters
    ----------
    pclDIct : Dictionary
        keys = drug family name, values = drugs belonging to that drug family
    metric : str
        query metric (eg, 'spearman' or 'wtcs')
    out : str
        output path
    summly_out_prefix : str
        directory prefix for summly output   
    pairwise_prefix : str
        directory prefix for graph output
    cell_match_mode : bool
        if True - run summly using cell line specific mode
        if False - run summly in non-matched mode
    row_space : str
        probe row_space
    '''
    def __init__(self, pclDict, 
                        metric, 
                        out,
                        pairwise_prefix='pairwise_matrices',
                        rankpt_metric='mean_rankpt_4',
                        sum_score_metric='sum_score_4',
                        row_space = 'lm',
                        cell_match_mode=True):
        '''
        Initialize a new instance of pathway connection analysis
        
        '''
        # set output directories
        self.out = out
        if not path.exists(self.out):
            os.mkdir(self.out)
        # set dir    
        self.pairwise_mtrx_out = out + '/' + pairwise_prefix
        if not path.exists(self.pairwise_mtrx_out):
            os.mkdir(self.pairwise_mtrx_out)
        # set query parameters
        self.pclDict = pclDict
        self.metric = metric
        self.row_space = row_space
        self.cell_match_mode = cell_match_mode
        self.sum_score_metric= sum_score_metric
        self.rankpt_metric = rankpt_metric
        # make set of all input compounds
        cpList = [item for sublist in pclDict.values() for item in sublist]
        cpSet = set(cpList)
        self.cpSet = cpSet

    def get_inames(self):
        '''
        get pert_inames for each input compound
        '''
        cm = mu.CMapMongo()
        inameID = cm.find({'pert_id':{'$in':list(self.cpSet)}},
                    {'pert_id':True,'pert_iname':True},
                    toDataFrame=True)
        inameSer = pd.Series(data=inameID['pert_iname'])
        inameSer.index = inameID['pert_id']
        inameDict = inameSer.to_dict()
        self.inameDict = inameDict

    def load_summly_mtrx(self):
        '''
        load pre-calculated summly matrix

        Parameters
        ----------

        '''
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
        rnkptNoDup = rnkptNoDup.reindex(index=colShortSet,columns=colShortSet)
        self.rnkpt_matrix_full = rnkptNoDup
        dictVals = self.pclDict.values()
        dictArray = [item for sublist in dictVals for item in sublist]
        dictUnique = list(set(dictArray))
        cpTested = [x for x in dictUnique if x in self.rnkpt_matrix_full.index.values]
        self.cpTested = cpTested
        testedMtrx = self.rnkpt_matrix_full.reindex(index=cpTested,columns=cpTested)
        self.rnkpt_matrix = testedMtrx
        ### dictionary of tested results
        pclResultDict = {}
        for grpName in self.pclDict:
            grpRes = [cp for cp in self.pclDict[grpName] if cp in cpTested] # leave out compounds that don't have summly data
            if len(grpRes) > 1:
                pclResultDict[grpName] = grpRes
        self.pclResultDict = pclResultDict
        ### make percent_summly matrix
        # cp_percent_summlyMtrx = pd.DataFrame()
        # for icp in iColFirst:
        #     cpSer = rnkpt.ix[:,icp]
        #     cpSer = cpSer[~np.isnan(cpSer)]
        #     cpSer = cpSer.order(ascending=False)
        #     nTested = len(cpSer)
        #     cpFrmRank = pd.DataFrame((np.arange(nTested)+1),
        #                 index=cpSer.index,
        #                 columns=['rank'])
        #     cpPerc = cpFrmRank/(float(nTested)+1) # percent rank of summly results
        #     pd.concat([cp_percent_summlyMtrx,cpPerc],
        #             ignore_index=False,
        #             axis=1)

    def test_groups(self,
        make_heatmaps=True,
        group_size_min=3):
        '''
        test the inter-connection of comound classes - make plots
        
        '''
        sumScoreDict = {} #matrix of connections for each group
        percSummlyDict = {}
        rnkptDict = {}
        for grpName in self.pclResultDict:
            if grpName == '-666':
                continue
            grp = self.pclResultDict[grpName]
            if not grp: # skip if grp is empty
                continue
            if len(grp) < group_size_min:
                continue
            grpInames = [self.inameDict[cp] for cp in grp] # inames
            grpZip = zip(*[grp,grpInames])
            # matrices for group connections
            if self.cell_match_mode:
                [grp_sum_score, 
                grp_rankpt, 
                grp_PercSummly] = self.pairwise_calc_match(grp)
            else:
                [grp_sum_score, 
                grp_rankpt, 
                grp_PercSummly] = self.pairwise_calc_nonmatch(grp)
            ### take averages of the upper and lower matrix segments
            av_grp_sum_score = self.av_mtrx(grp_sum_score)
            av_grp_rankpt = self.av_mtrx(grp_rankpt)
            av_grp_PercSummly = self.av_mtrx(grp_PercSummly)
            # av_grp_rank = av_mtrx(grp_rank)
            ### skip if rankpoint matrix is empty
            boolNan = np.isnan(av_grp_rankpt) 
            if boolNan.all(): #if all ranpoint matrix is nan skip
                continue            
            ### store averages (upper only)
            sumScoreDict[grpName] = self.av_mtrx_upper(grp_sum_score)
            rnkptDict[grpName] = self.av_mtrx_upper(grp_rankpt)
            percSummlyDict[grpName] = self.av_mtrx_upper(grp_PercSummly)
            ### write matrices
            rpF = os.path.join(self.pairwise_mtrx_out,grpName + '_mean_rankpt_matrix')
            self.write_pairwise_mtrx(grpZip,av_grp_rankpt,rpF)
            ### write rankpoint matrix
            psF = os.path.join(self.pairwise_mtrx_out,grpName + '_average_percent_summly_matrix')
            self.write_pairwise_mtrx(grpZip,av_grp_PercSummly,psF)
            ### plotting
            if make_heatmaps:
                outF = os.path.join(self.pairwise_mtrx_out,grpName + '_compound_group_heatmap.png')
                self.make_group_heatmap(grpName, grp, av_grp_rankpt, av_grp_PercSummly, outF)
        self.sumScoreDict = sumScoreDict
        self.rnkptDict = rnkptDict
        self.percSummlyDict = percSummlyDict

    def check_shRNA_connection(self):
        '''
        -if the genetic target of a compound has been knocked down, 
        check the connection between drug and gene KD
        -keys to input pclDict must be gene symbols

        Parameters
        ----------

        '''
        # make dictionary of CGS pert_ids to gene symbols
        CM = mu.CMapMongo()
        CGSpertIDs = CM.find({'pert_type':'trt_sh.cgs'},
                    {'pert_iname':True,'pert_id':True},
                    toDataFrame=True)
        CGSpertIDs.index = CGSpertIDs['pert_id']
        CGSser = CGSpertIDs['pert_iname']
        CGSpertDict = CGSser.to_dict()
        # check whick targets have KD and are in summly space
        summSpace = self.rnkpt_matrix_full.index
        pertIDs = CGSpertDict.keys()
        inames = CGSpertDict.values()
        summCGS = [x for x in summSpace if CGSpertDict.has_key(x)]
        cgsSer = pd.Series(CGSpertDict)
        # Series of pert_ids / inames in summly space
        cgsSer = cgsSer.reindex(summCGS)
        pert_to_iname = cgsSer.to_dict()
        iname_to_pert = {v:k for k, v in pert_to_iname.items()}
        # which targets have KDs in summspace
        targetList = self.pclResultDict.keys()
        targetsWithKD = [x for x in targetList if x in cgsSer.values]
        rnkptDict = {}
        rnkptExtendedList = []
        for tarGene in targetsWithKD:
            cpsTargetd = self.pclResultDict[tarGene]
            genePert = iname_to_pert[tarGene]
            rnkptList = {}
            for cp in cpsTargetd:
                rnkptVal = self.rnkpt_matrix_full.ix[genePert,cp]
                if not np.isnan(rnkptVal):
                    rnkptList[cp] = rnkptVal
                    rnkptExtendedList.append(rnkptVal)
            rnkptDict[tarGene] = rnkptList
        plt.hist(rnkptExtendedList,30)
        plt.xlim([-100,100])
        plt.ylabel('freq',fontweight='bold')
        plt.xlabel('rank point',fontweight='bold')
        plt.title('connection of compounds to their expected molecular target')
        outF = os.path.join(self.pairwise_mtrx_out,'target_id_hist.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        #write connection results to a file
        targetIDFrm = pd.DataFrame()
        for gene in rnkptDict:
            dGene = rnkptDict[gene]
            dSer = pd.Series(dGene)
            dSer.name = 'rnkpt_val'
            dFrm = pd.DataFrame(dSer)
            dFrm['target_gene'] = gene
            targetIDFrm = pd.concat([targetIDFrm,dFrm],axis=0)
        inameSer = pd.Series(self.inameDict)
        inameMtch = inameSer.reindex(targetIDFrm.index)
        targetIDFrm['iname'] = inameMtch.values
        tIDfile = self.pairwise_mtrx_out + '/target_id_summary.txt'
        targetIDFrm.to_csv(tIDfile)

    # def make_connection_hist(self,path):
    #     '''
    #     make histograms to illustrate how class connections apear 
    #     on the top of the summly rank
    #     '''
    #     for grpName in self.pclResultDict:
    #         grp = self.pclResultDict[grpName]
    #         #make symetric matrix of 
    #         for brd in grp:
    #             # find brd_tp in columns
    #             for brd_tp in self.brdTpDict[brd]:
    #                 resSer = self.rankptInFrm[brd_tp]
    #                 resSer.ix[grp]
    #             #look at each column
                

    #     rankpt_metric='mean_rankpt_4'
    #     summDir = [f for f in os.listdir(path) if os.path.isdir(path+'/'+f)]
    #     # resBRD = [x[:13] for x in summDir]
    #     # resBRDset = set(resBRD)
    #     summPath = [path+'/'+f for f in summDir]

    def get_sig_ids(self,write_grps=True):
        '''
        identify which compounds which are in summly space - for each compound
        write sig_ids to a file
        '''
        sumSpaceFile = '/xchip/cogs/sig_tools/sig_query/results/a2_internal_wtcs.lm/query_info_n73597.txt'
        summSpace = pd.read_csv(sumSpaceFile,sep='\t')
        notSummSpace = []
        isSummSpace = []
        sigIDdict = {}
        for brd in self.cpSet:
            if not brd in summSpace['pert_id'].values:
                notSummSpace.append(brd) # not enough sigatures for summly
            else:
                isSummSpace.append(brd)
                brdSpace = summSpace[summSpace['pert_id'] == brd]
                brdGrouped = brdSpace.groupby('cell_id')
                brdSigs = brdGrouped['sig_id'].first()
                sigIDdict[brd] = brdSigs.values
                if write_grps:
                    brdFile = self.summly_out + '/' + brd + '.grp'
                    brdSigs.to_csv(brdFile,index=False)
                ### if multiple signatures exist for a signature in a cell line
                ### choose the one with the highest replicate correlation            
                # ccGrouped = brdGrouped['distil_cc_q75']
                # groups = brdGrouped.groups
                # #loop through each cell line
                # for cell in groups:
                #     icell = groups[cell]
                #     g = brdSpace.ix[icell]['distil_cc_q75']
                # brdSingles 
        self.notSummSpace = notSummSpace
        self.isSummSpace = isSummSpace
        self.sigIDdict = sigIDdict

    def run_summly(self,rerun_mode=False):
        '''
        run summly on brd files generated 
        
        Parameters
        ----------
        rerun_mode : bool
            if True - re-run those perts who are in the job falure list from a previous run
            if False - run all input perts in summly space
        '''
        # compounds to run
        if rerun_mode:
            submitList = self.noSummlyResult
        else:
            submitList = self.isSummSpace
        # generate and run summly command
        for brd in submitList:
            summlyMtrx = '/xchip/cogs/sig_tools/sig_query/results/a2_internal_wtcs.lm/'
            querySpace = self.summly_out + '/' + brd + '.grp'
            cmd = ' '.join(['rum -q hour',
                     '-d sulfur_io=100',
                     '-o ' + self.summly_out,
                     '-x sig_summly_tool ' + summlyMtrx,
                     '--query_space ' + querySpace,
                     '--group_query ' + self.cell_match_mode_str,
                     '--out ' + self.summly_out])
            os.system(cmd)
       
    def make_summly_path_dict(self,path):
        '''
        take a path to summly results - make a dictionary of paths for every compound
        '''
        summDir = [f for f in os.listdir(path) if os.path.isdir(path+'/'+f)]
        summPath = [path+'/'+f for f in summDir]
        #put the path for each cp in a dictionary
        cpPathDict = {}
        jobFalureList = []
        for path in summPath:
            cpDirs = [f for f in os.listdir(path) if os.path.isdir(path+'/'+f)]
            if not cpDirs:
                continue
            if self.cell_match_mode:
                cp = cpDirs[0]
                # make sure this cp is one being tested
                if cp in self.cpSet:
                    cpPathDict[cp] = path + '/' + cp
            else:
                brdFull = cpDirs[0].split('_')[3]
                cp = brdFull[:13]
                sigDict = {}
                for sig in cpDirs:
                    sigDict[sig] = path + '/' + sig
                # make sure this cp is one being tested
                if cp in self.cpSet:
                    cpPathDict[cp] = sigDict
        #check for job failures - clear out failed directories
        if self.cell_match_mode:
            for cp in cpPathDict:
                sumFile = cpPathDict[cp] + '/'+cp + '_summly.txt'
                if not os.path.isfile(sumFile):
                    jobFalureList.append(cp)
                    #remove direcotry for failed jobs
                    cpDir = cpPathDict[cp] 
                    sDir = os.path.dirname(cpDir) #job direcotry
                    shutil.rmtree(sDir)
        else:
            for cp in cpPathDict:
                sigs = self.sigIDdict[cp]
                for sig in sigs:
                    jobIncomplete = False
                    if not cpPathDict[cp].has_key(sig):
                        jobIncomplete = True
                    else:
                        path1 = cpPathDict[cp][sig]
                        sumFile = path1 + '/'+sig + '_summly.txt'
                        if not os.path.isfile(sumFile):
                            jobIncomplete = True
                    if jobIncomplete:
                        jobFalureList.append(cp)
                        #remove direcotry for failed jobs
                        if cpPathDict.has_key(cp):
                            if cpPathDict[cp].has_key(sig):
                                cpDir = cpPathDict[cp][sig]
                                sDir = os.path.dirname(cpDir) #job direcotry
                                try:
                                    print 'incomplete deleting: ' + sDir
                                    shutil.rmtree(sDir)
                                except OSError:
                                    continue
        # identify which cps have summly result
        summlyLeftOut = []
        for brd in self.isSummSpace:
            if not brd in cpPathDict:
                summlyLeftOut.append(brd)
        #new pclDict with only summSpace members
        pclResultDict = {}
        for grpName in self.pclDict:
            grpRes = [cp for cp in self.pclDict[grpName] if cp in cpPathDict] # leave out compounds that don't have summly data
            pclResultDict[grpName] = grpRes
        self.pclResultDict = pclResultDict
        self.cpPathDict = cpPathDict
        self.noSummlyResult = summlyLeftOut
        self.jobFalureList = jobFalureList

    # def check_incomplete_jobs(self,path):
    #     '''
    #     check for incomplete jobs
    #     '''
    #     for cp in self.isSummSpace:
    #         sigs = self.sigIDdict[cp]
    #         for sig in sigs:
    #             jobIncomplete = False
    #             if not self.cpPathDict[cp].has_key(sig):
    #                 print cp + ' : ' + sig + ' - incomplete '
    #                 jobIncomplete = True
    #             else:
    #                 path1 = cpPathDict[cp][sig]
    #                 sumFile = path1 + '/'+sig + '_summly.txt'
    #                 if not os.path.isfile(sumFile):
    #                     jobIncomplete = True

    def test_DOS_against_groups(self,
        make_heatmaps=True,
        group_size_min=3,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4'):
        '''
        -test how DOS compounds relate to known PCL groups
        -Do symetric look ups for both query directions for a given 
        pairwise relationship  
        '''
        dos_rnkpt = {}
        for doscp in self.pclDict['DOS']:
            dos_rnkpt_dist = {}
            for grpName in self.pclDict: # top connecting cps?
                if (grpName == '-666') or (grpName == '-666'):
                    continue
                grp = self.pclDict[grpName]
                grp.append(doscp) # add dos cp to group and test for relatedness 
                if not grp: # skip if grp is empty
                    continue
                grp = [cp for cp in grp if cp in self.cpPathDict] # leave out compounds that don't have summly data
                if len(grp) < group_size_min:
                    continue
                grpInames = [self.inameDict[cp] for cp in grp] # inames
                grpZip = zip(*[grp,grpInames])
                # matrices for group connections
                if self.cell_match_mode:
                    [grp_sum_score, 
                    grp_rankpt, 
                    grp_PercSummly] = self.pairwise_calc_match(grp,
                                                    sum_score_metric,
                                                    rankpt_metric)
                else:
                    [grp_sum_score, 
                    grp_rankpt, 
                    grp_PercSummly] = self.pairwise_calc_nonmatch(grp,
                                                    sum_score_metric,
                                                    rankpt_metric)
                ### take averages of the upper and lower matrix segments
                av_grp_sum_score = self.av_mtrx(grp_sum_score)
                av_grp_rankpt = self.av_mtrx(grp_rankpt)
                av_grp_PercSummly = self.av_mtrx(grp_PercSummly)
                dos_rnkpt_dist[grpName] = av_grp_rankpt[:,-1]
            dos_rnkpt[doscp] = dos_rnkpt_dist
            # set up graphs
            rnkptSumList = []
            tickList = []
            for gName in dos_rnkpt_dist:
                # rankpt
                flatM2 = dos_rnkpt_dist[gName]
                flatM2 = flatM2[~np.isnan(flatM2)] # remove nan
                rnkptSumList.append(flatM2)
                #names
                tickList.append(gName)
            #sumscore boxplot
            plt.boxplot(rnkptSumList,vert=0)
            plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
            plt.tick_params(labelsize=8)
            plt.ylabel('compound class',fontweight='bold')
            plt.xlabel('rank point',fontweight='bold')
            plt.title('DOS cp: ' + doscp + ' - distribution of rank point values by group',fontweight='bold')
            outF = os.path.join(self.pairwise_mtrx_out,doscp+ '_DOS_rankpt_boxplot.png')
            plt.savefig(outF, bbox_inches='tight',dpi=200)
            plt.close()
        self.dos_rnkpt_by_goup = dos_rnkpt

    def test_DOS_queries(self,
        make_heatmaps=True,
        make_boxplots=True, 
        rnkpt_thresh=90,       
        group_size_min=3,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4'):
        '''
        -test how DOS compounds relate to known PCL groups
        -only use query results for DOS compounds - non symetric lookup

        Parameters
        ----------
        rnkpt_thresh : int
            threshold for making heatmap of a particular dos-group connection - median rankpoint 

        '''
        dos_rnkpt = {}
        for doscp in self.pclDict['DOS']:
            dfQuery = pd.DataFrame()
            if self.cell_match_mode:
                inFile = '/'.join([self.cpPathDict[doscp],
                                doscp+'_summly.txt'])
                sumRes = pd.io.parsers.read_csv(inFile,sep='\t')
                sumRes = sumRes.replace(-666,np.nan)
                cpRes = sumRes[sumRes['pert_type'] == 'trt_cp']
                cpRes['rank'] = np.arange(1,len(cpRes)+1)
                cpRes.index = cpRes['pert_id']
                for grpName in self.pclDict: # top connecting cps?
                    if (grpName == '-666') or (grpName == '-DOS'):
                        continue
                    grp = self.pclDict[grpName]
                    if not grp: # skip if grp is empty
                        continue
                    grp = [cp for cp in grp if cp in self.cpPathDict] # leave out compounds that don't have summly data
                    if len(grp) < group_size_min:
                        continue
                    # load in dos group
                    grpRes = pd.DataFrame() #clear frame
                    grpRes = cpRes.ix[grp]
                    grpRes['DOS_query'] = doscp
                    grpRes['group_name'] = grpName
                    grpRes = grpRes.rename(columns={'pert_id':'queried_brd'})
                    dfQuery = pd.concat([dfQuery,grpRes],axis=0,ignore_index=True)
                    grpRnkpt = grpRes[rankpt_metric]
                    dos_rnkpt_dist[grpName] = grpRnkpt
                    ### make heatmap for best connectors
                    if make_heatmaps and (grpRnkpt.median() >rnkpt_thresh):
                        print doscp + ' strong connection with ' + grpName + ' group'
                        grpPlusDos = grp + [doscp]
                        grpInames = [self.inameDict[cp] for cp in grpPlusDos] # inames
                        grpZip = zip(*[grpPlusDos,grpInames])
                        # matrices for group connections
                        if self.cell_match_mode:
                            [grp_sum_score, 
                            grp_rankpt, 
                            grp_PercSummly] = self.pairwise_calc_match(grpPlusDos,
                                                            sum_score_metric,
                                                            rankpt_metric)
                        else:
                            [grp_sum_score, 
                            grp_rankpt, 
                            grp_PercSummly] = self.pairwise_calc_nonmatch(grpPlusDos,
                                                            sum_score_metric,
                                                            rankpt_metric)
                        ### take averages of the upper and lower matrix segments
                        av_grp_sum_score = self.av_mtrx(grp_sum_score)
                        av_grp_rankpt = self.av_mtrx(grp_rankpt)
                        av_grp_PercSummly = self.av_mtrx(grp_PercSummly)
                        ### plotting
                        if make_heatmaps:
                            outF = os.path.join(self.pairwise_mtrx_out, doscp + '_' + grpName + '_comparison_heatmap.png')
                            self.make_group_heatmap(grpName, grpPlusDos, av_grp_rankpt, av_grp_PercSummly, outF)
            dos_rnkpt[doscp] = dfQuery
            if make_boxplots:
                grouped = dfQuery.groupby('group_name')
                dosMedians = grouped[rankpt_metric].median() # meadian ranpt for dos-group connections
                dosMedians.sort()
                dosMedians = dosMedians[dosMedians.notnull()]
                rnkptSumList = []
                tickList = []
                for gName in dosMedians.index:
                    # rankpt
                    flatM2 = dos_rnkpt_dist[gName]
                    flatM2 = flatM2[~np.isnan(flatM2)] # remove nan
                    rnkptSumList.append(flatM2)
                    #names
                    tickList.append(gName)
                #sumscore boxplot
                plt.boxplot(rnkptSumList,vert=0)
                plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
                plt.tick_params(labelsize=8)
                plt.ylabel('compound class',fontweight='bold')
                plt.xlabel('rank point',fontweight='bold')
                plt.xlim((-100,100))
                plt.title('DOS cp: ' + doscp + ' - distribution of rank point values by group',fontweight='bold')
                outF = os.path.join(self.pairwise_mtrx_out,doscp+ '_DOS_rankpt_boxplot.png')
                plt.savefig(outF, bbox_inches='tight',dpi=200)
                plt.close()
        self.dos_rnkpt_by_goup = dos_rnkpt

    def test_class_interrelatedness_old(self,
        make_heatmaps=True,
        make_boxplots=True,
        make_group_by_cp_mtrx=True,
        rankpt_metric='mean_rankpt_4'):
        '''
        -test the interrelated each of the PCL classes 

        Parameters
        ----------
        rnkpt_metric : str
            summly metric on which to sort
        make_heatmap : bool
            make summary graph
        make_boxplot : bool
            make boxplot for each groups
        make_group_by_cp_mtrx : bool
            matrix of median summly results - group x summly compounds

        '''
        dfQuery = pd.DataFrame()
        for brd in self.isSummSpace:
            grp = self.isSummSpace[:]
            grp.remove(brd)
            if self.cell_match_mode:
                inFile = '/'.join([self.cpPathDict[brd],
                                brd+'_summly.txt'])
                sumRes = pd.io.parsers.read_csv(inFile,sep='\t')
                sumRes = sumRes.replace(-666,np.nan)
                cpRes = sumRes[sumRes['pert_type'] == 'trt_cp']
                cpRes['rank'] = np.arange(1,len(cpRes)+1)
                cpRes.index = cpRes['pert_id']
                # where does every other connection sit on the rank list?
                grpRes = pd.DataFrame() #clear frame
                grpRes = cpRes.ix[grp]
                grpRes['query_cp'] = brd
                # grpRes['group_name'] = grpName
                grpRes = grpRes.rename(columns={'pert_id':'queried_brd'})
                dfQuery = pd.concat([dfQuery,grpRes],axis=0,ignore_index=True)
        dfQuery.index = dfQuery['queried_brd']
        self.dfQuery = dfQuery
        ### make empty group-group relatedness frame
        interGrpFrm = pd.DataFrame()
        grpBycpFrm = pd.DataFrame()
        for grp1Name in self.pclResultDict:
            grp1Frm = pd.DataFrame()
            # for each group make a frame with connection data to all other groups
            grp1 = self.pclResultDict[grp1Name]
            grp1Set = set(grp1)
            for brd1 in grp1:
                brd1Conn = dfQuery[dfQuery['query_cp'] == brd1]
                for grp2Name in self.pclResultDict:
                    grp2Tmp = set(self.pclResultDict[grp2Name])
                    #remove compound from group2 if it is in group1
                    grp2 = grp2Tmp.difference(grp1Set)
                    if grp2:
                        grpRes = brd1Conn.ix[grp2]
                        grpRes['query_group'] = grp1Name
                        grpRes['queried_group'] = grp2Name
                        grp1Frm = pd.concat([grp1Frm,grpRes],axis=0,ignore_index=True)
            # organize by group
            grouped = grp1Frm.groupby('queried_group')
            grpMedians = grouped[rankpt_metric].median() # meadian ranpt for dos-group connections
            grpMedians.name = grp1Name
            medFrame = pd.DataFrame(grpMedians)
            interGrpFrm = pd.concat([interGrpFrm,medFrame],
                                axis=1,
                                ignore_index=False)
            # organize by compound
            cpGrouped = grp1Frm.groupby('queried_brd')
            cpGrpMedians = cpGrouped[rankpt_metric].median() # meadian ranpt for dos-group connections
            cpGrpMedians.name = grp1Name
            cpMedFrame = pd.DataFrame(cpGrpMedians)
            grpBycpFrm = pd.concat([grpBycpFrm,cpMedFrame],
                                axis=1,
                                ignore_index=False)
            if make_boxplots:
                grpMedians.sort()
                grpMedians = grpMedians[grpMedians.notnull()]
                rnkptSumList = []
                tickList = []
                for gName in grpMedians.index:
                    # rankpt
                    iGrp = grouped.groups[gName]
                    grpConn = grp1Frm.ix[iGrp]
                    flatM2 = grpConn[rankpt_metric].values # all pairwise values
                    querySlice = grpConn.groupby('query_cp')
                    # median pairwise relationship for each compound in grp1
                    flatM2 = querySlice.median()[rankpt_metric].values
                    flatM2 = flatM2[~np.isnan(flatM2)] # remove nan
                    rnkptSumList.append(flatM2)
                    #names
                    tickList.append(gName)
                #sumscore boxplot
                plt.boxplot(rnkptSumList,vert=0)
                plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
                plt.tick_params(labelsize=8)
                plt.ylabel('compound class',fontweight='bold')
                plt.xlabel('rank point',fontweight='bold')
                plt.xlim((-100,100))
                plt.title(grp1Name + ' - distribution of inter-group rank point values',fontweight='bold')
                outF = os.path.join(self.pairwise_mtrx_out,grp1Name+ '_inter_group_rankpt_boxplot.png')
                plt.savefig(outF, bbox_inches='tight',dpi=200)
                plt.close()
        # cluster inter-groups Frame and make heatmap
        self.interGrpFrm = interGrpFrm
        self.grpBycpFrm = grpBycpFrm
        if make_heatmaps:
            self.inter_group_cluster(interGrpFrm)
        # make group x compound matrix
        if make_group_by_cp_mtrx:
            mtrxFile = self.pairwise_mtrx_out + '/group_by_compound_rankpt_mtrx.txt'
            grpBycpFrm.to_csv(mtrxFile)

    def test_class_interrelatedness(self,
        make_heatmaps=True,
        make_boxplots=True,
        make_group_by_cp_mtrx=True):
        '''
        -test the interrelated each of the PCL classes 

        Parameters
        ----------
        make_heatmap : bool
            make summary graph
        make_boxplot : bool
            make boxplot for each groups
        make_group_by_cp_mtrx : bool
            matrix of median summly results - group x summly compounds

        '''
        ### make empty group-group relatedness frame
        interGrpFrm = pd.DataFrame()
        grpBycpFrm = pd.DataFrame()
        for grp1Name in self.rankpt_group_medians.index:
            grp1Frm = pd.DataFrame()
            grpMedianDict = {}
            # for each group make a frame with connection data to all other groups
            grp1 = self.pclResultDict[grp1Name]
            # grp1arbitraryTps = [self.brdTpDict[brd][0] for brd in grp1]
            grp1Set = set(grp1)
            # for brd1 in grp1:
            for grp2Name in self.rankpt_group_medians.index:
                grp2Tmp = set(self.pclResultDict[grp2Name])
                #remove compound from group2 if it is in group1
                grp2 = grp2Tmp.difference(grp1Set)
                if grp2:
                    grpRes = self.rnkpt_matrix.ix[grp2,grp1]
                    medianRnkpt = stats.nanmean(grpRes.unstack())
                    grpRes['query_group'] = grp1Name
                    grpRes['queried_group'] = grp2Name
                    grp1Frm = pd.concat([grp1Frm,grpRes],axis=0,ignore_index=True)
                    grpMedianDict[grp2Name] = medianRnkpt
            # organize by group
            grouped = grp1Frm.groupby('queried_group')
            # # grpMedians = grouped[rankpt_metric].median() # meadian ranpt for dos-group connections
            # grpMedians = grouped.median()
            grpMedians = pd.Series(grpMedianDict)
            grpMedians.name = grp1Name
            medFrame = pd.DataFrame(grpMedians)
            interGrpFrm = pd.concat([interGrpFrm,medFrame],
                                axis=1,
                                ignore_index=False)
            # organize by compound
            # cpGrouped = grp1Frm.groupby('queried_brd')
            # cpGrpMedians = cpGrouped[rankpt_metric].median() # meadian ranpt for dos-group connections
            # cpGrpMedians.name = grp1Name
            # cpMedFrame = pd.DataFrame(cpGrpMedians)
            # grpBycpFrm = pd.concat([grpBycpFrm,cpMedFrame],
            #                     axis=1,
            #                     ignore_index=False)
            if make_boxplots:
                grpMedians.sort()
                grpMedians = grpMedians[grpMedians.notnull()]
                rnkptSumList = []
                tickList = []
                for gName in grpMedians.index:
                    # rankpt
                    iGrp = grouped.groups[gName]
                    # grpConn = grp1Frm.ix[iGrp]
                    flatM2 = grp1Frm.ix[iGrp,grp1].values
                    flatM2 = flatM2[~np.isnan(flatM2)] # remove nan
                    # flatM2 = grpConn[rankpt_metric].values # all pairwise values
                    # querySlice = grpConn.groupby('query_cp')
                    # # median pairwise relationship for each compound in grp1
                    # flatM2 = querySlice.median()[rankpt_metric].values
                    # flatM2 = flatM2[~np.isnan(flatM2)] # remove nan
                    rnkptSumList.append(flatM2)
                    #names
                    tickList.append(gName)
                #sumscore boxplot
                plt.boxplot(rnkptSumList,vert=0)
                plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
                plt.tick_params(labelsize=8)
                plt.ylabel('compound class',fontweight='bold')
                plt.xlabel('rank point',fontweight='bold')
                plt.xlim((-100,100))
                plt.title(grp1Name + ' - distribution of inter-group rank point values',fontweight='bold')
                outF = os.path.join(self.pairwise_mtrx_out,grp1Name+ '_inter_group_rankpt_boxplot.png')
                plt.savefig(outF, bbox_inches='tight',dpi=200)
                plt.close()
        # cluster inter-groups Frame and make heatmap
        interGrpFrm = interGrpFrm.reindex(index=self.rankpt_group_medians.index, 
                                        columns=self.rankpt_group_medians.index)
        self.interGrpFrm = interGrpFrm
        # self.grpBycpFrm = grpBycpFrm
        if make_heatmaps:
            self.inter_group_cluster(interGrpFrm)
        # make group x compound matrix
        # if make_group_by_cp_mtrx:
        #     mtrxFile = self.pairwise_mtrx_out + '/group_by_compound_rankpt_mtrx.txt'
        #     grpBycpFrm.to_csv(mtrxFile)
        # make histogram of inter-group median connections
        interMed = interGrpFrm.unstack().values
        interMed = interMed[~np.isnan(interMed)]
        # 'median ranpoint connection between groups'

    def inter_group_cluster(self,
        interGrpFrm):
        '''
        -Takes a matrix of median pairwise summly results - group x group
        -clusters this matrix and makes heatmap
        
        Parameters
        ----------
        interGrpFrm : Pandas DataFrame
            -a matrix of median summly results
            -n by n 
            -index and columns are the PCL names

        '''
        scoreFrm = pd.DataFrame(interGrpFrm.copy()/float(100),columns=interGrpFrm.index)
        meanScore = stats.nanmean(scoreFrm.values.ravel())
        np.fill_diagonal(scoreFrm.values, 1)
        scoreFrm = scoreFrm.fillna(meanScore)
        annots = pd.DataFrame(data=interGrpFrm.index.values,index=interGrpFrm.index,columns=['group_name'])
        hcl = HClust(scoreFrm.fillna(0), annots, out = self.pairwise_mtrx_out, 
                        symmfun = np.max, cutoff = 0.5, 
                        convert_to_distance=True)
        hcl.cluster()
        # hcl.draw_heatmap(vmin = -1, vmax = 1,
        #                   cmap = cm.RdBu_r,
        #                   col_label = 'group_name',
        #                   row_label = None,
        #                   showfig=True,
        #                   title = 'summly clustering')        
        clust_order = hcl.cluster_order.order().index
        annots['hclust_assignment'] = hcl.cluster_assignment
        rankptClustered = interGrpFrm.reindex(clust_order,columns=clust_order)
        meanedRnkptClust = self.av_mtrx(rankptClustered.values)
        outF = os.path.join(self.pairwise_mtrx_out,'median_rankpt_class_class_clustering.png')
        fig = plt.figure(1, figsize=(20, 8))
        colors.set_color_map()
        plt.suptitle('median PCL group interrelatedness',fontsize=14, fontweight='bold')
        plt.title('mean_rankpt_4')
        plt.imshow(meanedRnkptClust,
                interpolation='nearest',
                vmin=-100, 
                vmax=100)
                # cmap=matplotlib.cm.RdBu_r,        
        ytcks = [x for x in clust_order]
        plt.xticks(np.arange(len(rankptClustered))-1, ytcks,rotation=75)
        plt.yticks(np.arange(len(rankptClustered)),ytcks)
        plt.colorbar()
        fig.savefig(outF, bbox_inches='tight')
        plt.close()
        self.group_clust_order = clust_order

    def inter_group_line_graph(self,
        rankpt_thresh=50):
        '''
        -make a graph where the line thickness represents strength of group-group
        interaction
        
        Parameters
        ----------
        interGrpFrm : Pandas DataFrame
            -a matrix of median summly results
            -n by n 
            -index and columns are the PCL names
        interGrpFrm : list
            list to represent label orders

        '''
        # upper vs full matrix
        # symFrm = self.av_mtrx_upper(self.interGrpFrm.values)
        # symFrm = self.av_mtrx(self.interGrpFrm.values)
        # reorder matrix acording to cluster results
        clusteredMtrx = self.interGrpFrm.reindex(index=self.group_clust_order,
                                columns=self.group_clust_order)
        symFrm = self.av_mtrx(clusteredMtrx.values)
        num_groups = self.interGrpFrm.shape[0]
        lwList = []
        fig = plt.figure(1, figsize=(11, 11))
        # loop through each coordinate of the matrix
        for i in range(num_groups):
            for j in range(num_groups):
                pairVal = symFrm[i,j]
                if np.abs(pairVal) > rankpt_thresh:
                    # lw = np.power(pairVal/float(85),10)
                    # lw = np.power(pairVal/float(85),4)
                    lw = np.power(pairVal/float(65),4)
                    lwList.append(lw)
                    plt.plot([1,2],[i,j],'b',linewidth=lw)
                    # plt.plot([2,1],[j,i],'b',linewidth=lw)
        ytcks = [x.replace('-inhibitor','') for x in self.interGrpFrm.index]
        plt.xlim((1,2))
        plt.tick_params(labelright=True,bottom=False)
        plt.yticks(np.arange(len(ytcks)),ytcks)
        outF = self.pairwise_mtrx_out + '/inter_group_line_graph.png'
        fig.savefig(outF, bbox_inches='tight')
        plt.close()

    def group_probe_frq_plot(self,
        make_heatmaps=True,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4'):
        '''
        test relative occurance of up/dn regulation of probes for a specific group
        
        '''
        brd = 'BRD-K02130563'
        sigs = po.sigIDdict[brd]
        sig = sigs[0]
        #
        afPath = cmap.score_path
        gt = gct.GCT()
        gt.read(src=afPath,cid=sigs,rid='lm_epsilon')
        zFrm = gt.frame
        # zFrm = pd.DataFrame(data=gt.matrix,
        #                     index=gt.get_rids(),
        #                     columns=sigs)
        # take modz of signature group
        modZed = modzsig.modzsig(zFrm)
        modZed = modZed.order()
        #pick a group 
        # grpName = 'tubulin'
        grpName = 'HDAC-inhibitor'
        #get all sig_ids for that group
        grpSigList = []
        for brd in self.pclResultDict[grpName]:
            grpSigList.extend(self.sigIDdict[brd])
        #query for up/dn probes
        cm = mu.CMapMongo()
        regFrm = cm.find({'sig_id':{'$in':list(grpSigList)}},
                    {'sig_id':True,'pert_id':True,'pert_iname':True,'up50_lm':True,'dn50_lm':True},
                    toDataFrame=True)
        # count dn probe freq
        nInstances = regFrm.shape[0]
        dnNested = regFrm['dn50_lm'].values
        dnArray = [item for sublist in dnNested for item in sublist]
        dnSer = pd.Series(dnArray)
        dnCounts = dnSer.value_counts()
        zDnCounts = dnCounts.reindex_like(modZed)
        # count dn probe freq
        upNested = regFrm['up50_lm'].values
        upArray = [item for sublist in upNested for item in sublist]
        upSer = pd.Series(upArray)
        upCounts = upSer.value_counts()
        zUpCounts = upCounts.reindex_like(modZed)
        # adjust marker size
        upPercMkrs = np.divide(zUpCounts,nInstances) #divide by total instances to make for relative frequency
        dnPercMkrs = np.divide(zDnCounts,nInstances)
        upMkrs = np.multiply(upPercMkrs,100)
        dnMkrs = np.multiply(dnPercMkrs,100)
        upMkrs = upMkrs.replace(np.nan,0)
        dnMkrs = dnMkrs.replace(np.nan,0)
        # make plot 
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # ax.plot(s,s,'b')
        for j,sl in enumerate(modZed):
            ax.plot(j,1,'r.',markersize=upMkrs[j],alpha=.25)
            ax.plot(j,1,'b.',markersize=dnMkrs[j],alpha=.25)

    def cluster_all_cps(self,
        make_heatmaps=True,
        sum_score_metric='sum_score_4',
        rankpt_metric='mean_rankpt_4'):
        '''
        test the inter-connection of comound classes - make plots
        
        '''
        allbrds = self.cpPathDict.keys()
        grpInames = [self.inameDict[cp] for cp in allbrds] # inames
        grpZip = zip(*[allbrds,grpInames])
        # matrices for group connections
        if self.cell_match_mode:
            [grp_sum_score, 
            grp_rankpt, 
            grp_PercSummly] = self.pairwise_calc_match(allbrds,
                                            sum_score_metric,
                                            rankpt_metric)
        else:
            [grp_sum_score, 
            grp_rankpt, 
            grp_PercSummly] = self.pairwise_calc_nonmatch(allbrds,
                                            sum_score_metric,
                                            rankpt_metric)
        ### take averages of the upper and lower matrix segments
        av_grp_sum_score = self.av_mtrx(grp_sum_score)
        av_grp_rankpt = self.av_mtrx(grp_rankpt)
        av_grp_PercSummly = self.av_mtrx(grp_PercSummly)
        # av_grp_rank = av_mtrx(grp_rank)
        ### write matrices
        rpF = os.path.join(self.pairwise_mtrx_out,'whole_mean_rankpt_matrix')
        self.write_pairwise_mtrx(grpZip,av_grp_rankpt,rpF)
        ### write rankpoint matrix
        psF = os.path.join(self.pairwise_mtrx_out,'whole_percent_summly_matrix')
        self.write_pairwise_mtrx(grpZip,av_grp_PercSummly,psF)
        ### plotting
        if make_heatmaps:
            outF = os.path.join(self.pairwise_mtrx_out,'whole_compound_group_heatmap.png')
            self.make_group_heatmap('all_compounds', allbrds, av_grp_rankpt, av_grp_PercSummly, outF)
        self.whole_sum_score = av_grp_sum_score
        self.whole_rankpt = av_grp_rankpt
        self.whole_PercSummly = av_grp_PercSummly
        perform_HClust(allbrds, av_grp_rankpt,av_grp_PercSummly)

    def perform_HClust(self,grp,rankptMtrx,percSummlyMtrx):
        '''
        calculate pairwise matrices during summly cell-line-matched mode

        Parameters
        ----------
        grp : list of strings
            list of compound names belonging to a group with summly results
        mtrx : numpy matrix
            numpy matrix of values to be clustered

        ''' 
        rankptFrm = pd.DataFrame(rankptMtrx,index=grp,columns=grp)
        percSummlyFrm = pd.DataFrame(percSummlyMtrx,index=grp,columns=grp)
        scoreFrm = pd.DataFrame(rankptMtrx/float(100),index=grp,columns=grp)
        inames = [str(self.inameDict[brd]) for brd in grp]
        annots = pd.DataFrame(data=inames,index=grp,columns=['pert_iname'])
        hcl = HClust(scoreFrm, annots, out = self.pairwise_mtrx_out, 
                        symmfun = np.max, cutoff = 0.5, 
                        convert_to_distance=True)
        hcl.cluster()
        hcl.get_cluster_order()
        # hcl.draw_dendrogram(orientation = 'right', showfig=False)
        # hcl.save_dendrogram()
        clust_order = hcl.cluster_order.order().index
        annots['hclust_assignment'] = hcl.cluster_assignment
        # hcl.draw_heatmap(order = clust_order,
        #                   vmin = 0, vmax = 1,
        #                   cmap = cm.OrRd,
        #                   col_label = 'pert_iname',
        #                   row_label = None,
        #                   annots_top = [('hclust_assignment', False)],
        #                   showfig=True,
        #                   title = 'summly clustering')
        rankptClustered = rankptFrm.reindex(clust_order,columns=clust_order)
        percSummClustered = percSummlyFrm.reindex(clust_order,columns=clust_order)
        # plt.imshow(rankptClustered.values,
        #     interpolation='nearest',
        #     cmap=cm.Greens_r,
        #     vmin=-100,
        #     vmax=100)
        outF = os.path.join(self.pairwise_mtrx_out,'whole_compound_clustered.png')
        self.make_group_heatmap('all_compounds', 
                                cluster_order, 
                                rankptClustered.values, 
                                percSummClustered.values, 
                                outF)

    def pairwise_calc_match(self,grp):
        '''
        calculate pairwise matrices during summly matched mode

        Parameters
        ----------
        grp : list of strings
            list of compound names belonging to a group with summly results
        Returns
        ----------
        3 pairwise relatedness matrices

        ''' 
        # retrieve the pairwise relations among compounds in a specific group
        nGrp = len(grp)
        rnkptPairwiseFrm = self.rnkpt_matrix.ix[grp,grp]
        # ssPairwiseFrm = self.sum_scoreMtrx.ix[grp,grp] # sum_score
        # psPairwiseFrm = self.cp_percent_summlyMtrx.ix[grp,grp] # percent summly
        # set ss and ps matrices to nan for now
        ssPairwiseFrm = np.zeros((nGrp,nGrp))
        ssPairwiseFrm[:] = np.nan
        psPairwiseFrm = np.zeros((nGrp,nGrp))
        psPairwiseFrm[:] = np.nan
        return [ssPairwiseFrm, rnkptPairwiseFrm.values, psPairwiseFrm]
        # return [ssPairwiseFrm.values, rnkptPairwiseFrm.values, psPairwiseFrm.values]

    def write_pairwise_mtrx(self,inames_zip, mtrx,out):
        '''
        Write a matrix to file

        Parameters
        ----------
        inames_zip : list of tuples
            brds paired with inames
        mtrx : numpy.ndarray
            matrix of data
        out : str
            output path - no file extension       
        ''' 
        Hindex = pd.MultiIndex.from_tuples(inames_zip, names=['brd', 'iname'])
        sumScoreFrm = pd.DataFrame(mtrx,index=Hindex,columns=Hindex)
        sumScoreFrm.to_csv(out+'.txt',sep='\t')
        gc = gct.GCT()
        gc.build_from_DataFrame(sumScoreFrm)
        gc.write(out)

    def make_group_heatmap(self,group_name, grp,rankpt_mtrx, sumRank_mtrx, out):
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
        plt.subplot(121)
        plt.title('mean_rankpt_4')
        colors.set_color_map()
        plt.imshow(rankpt_mtrx,
                interpolation='nearest',
                vmin=-100, 
                vmax=100)
        ytcks = [self.inameDict[x] for x in grp]
        plt.xticks(np.arange(len(grp)), ytcks,rotation=75)
        plt.yticks(np.arange(len(grp)),ytcks)
        plt.colorbar()
        plt.subplot(122)
        plt.title('percent_summly')
        plt.imshow(sumRank_mtrx,
                interpolation='nearest',
                cmap=matplotlib.cm.Greens_r,
                vmin=0,
                vmax=1)
        plt.xticks(np.arange(len(grp)), ytcks,rotation=75)
        plt.yticks(np.arange(len(grp)),ytcks)
        plt.colorbar()
        fig.savefig(out, bbox_inches='tight')
        plt.close()

    def make_summary_boxplot(self):
        '''
        plot the distribution of pairwise connections across
        all groups
        
        '''        
        ##order goups acording to mean sum_score
        # sumMedianDict = {}
        # for gName in self.sumScoreDict:
        #     sumMtrx = self.sumScoreDict[gName] # sum score matrix for group
        #     medianSum = np.median(sumMtrx[~np.isnan(sumMtrx)])
        #     sumMedianDict[gName] = medianSum
        # sumMedianSer = pd.Series(sumMedianDict)
        # sumMedianSer = sumMedianSer.order()
        #order goups acording to median rankpt value
        rpMedianDict = {}
        for gName in self.rnkptDict:
            rpMtrx = self.rnkptDict[gName] # sum score matrix for group
            medianRp = np.median(rpMtrx[~np.isnan(rpMtrx)])
            rpMedianDict[gName] = medianRp
        rpMedianSer = pd.Series(rpMedianDict)
        rpMedianSer = rpMedianSer.order()
        #make boxplot of all connections
        sumScoreList = []
        rnkptSumList = []
        percSummList = []
        tickList = []
        for gName in rpMedianSer.index:
            # sum score setup
            m1 = self.sumScoreDict[gName]
            flatM = m1.flatten()
            flatM = flatM[~np.isnan(flatM)] # remove nan
            sumScoreList.append(flatM)
            # rankpt
            m2 = self.rnkptDict[gName]
            flatM2 = m2.flatten()
            flatM2 = flatM2[~np.isnan(flatM2)] # remove nan
            rnkptSumList.append(flatM2)
            # percent summly setup
            m3 = self.percSummlyDict[gName]
            flatM3 = m3.flatten()
            flatM3 = flatM3[~np.isnan(flatM3)] # remove nan
            percSummList.append(flatM3)
            #names
            tickList.append(gName)
        #sumscore boxplot
        plt.boxplot(sumScoreList,vert=0)
        plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
        plt.tick_params(labelsize=8)
        plt.ylabel('compound class',fontweight='bold')
        plt.xlabel('sum_score_4',fontweight='bold')
        plt.title('distribution of sum scores values by group',fontweight='bold')
        outF = os.path.join(self.pairwise_mtrx_out,'sumScore_boxplot.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        plt.close()
        #percent summly boxplot 
        plt.boxplot(rnkptSumList,vert=0)
        plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
        plt.tick_params(labelsize=8)
        plt.ylabel('compound class',fontweight='bold')
        plt.xlabel('rank point',fontweight='bold')
        plt.title('distribution of rankpt values by group',fontweight='bold')
        outF = os.path.join(self.pairwise_mtrx_out,'rankpt_boxplot.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        plt.close()
        #percent summly boxplot 
        plt.boxplot(percSummList,vert=0)
        plt.yticks(np.arange(1,len(tickList)+1),tickList,rotation=0)
        plt.tick_params(labelsize=8)
        plt.ylabel('compound class',fontweight='bold')
        plt.xlabel('percent summly',fontweight='bold')
        plt.title('distribution of percent summly rank scores by group',fontweight='bold')
        outF = os.path.join(self.pairwise_mtrx_out,'percent_summly_rank_boxplot.png')
        plt.savefig(outF, bbox_inches='tight',dpi=200)
        plt.close()
        # self.summly_group_medians = sumMedianSer
        self.rankpt_group_medians = rpMedianSer

    def av_mtrx_upper(self,mtrx):
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

    def av_mtrx(self,mtrx):
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
