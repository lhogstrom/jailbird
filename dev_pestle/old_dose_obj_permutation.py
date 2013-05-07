#! /usr/bin/env python
'''
create a class to handel common info for dose analysis
'''

import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.progress as progress
import cmap.analytics.FDR as FDR
import scipy.stats as stats
import random

    def permutation_template(self,n_permutation=10000):
        '''
        creates null distribution of correlations 
        returns correlations, p-values, and FDR correction 

        must add_from_gct() before running
        '''
        prog = progress.DeterminateProgressBar('template matching')
        gcto = self.gct

        doses = gcto.get_column_meta('pert_dose')
        perts = gcto.get_column_meta('pert_id')
        rids = gcto.get_rids()

        examineList = self.perts_at_dose
        num_perts = len(examineList)
        # templateMatchInd = {} #nested dict for index of significant probes
        # pvecDictParametric = {}
        pvecDictEmpirical = {}
        corr_null_distribution = {}
        # probe_template_corrs = {}
        for icmpd,unique_pert in enumerate(examineList):
            prog.update('template match {0}'.format(unique_pert),icmpd,num_perts)
            cid_inds = [i for i,x in enumerate(perts) if unique_pert in x]
            pert_doses = [float(doses[x]) for x in cid_inds]
            tmp_tup = zip(pert_doses,cid_inds)
            tmp_tup.sort()
            pert_doses,cid_inds = zip(*tmp_tup)
            pert_data = gcto.matrix[:,cid_inds]
            template_names = ['linear', 'log10', 'log2']
            # templateMatchInd[unique_pert] = {}
            # pvecDictParametric[unique_pert] = {}
            pvecDictEmpirical[unique_pert] = {}
            corr_null_distribution[unique_pert] = {}
            # probe_template_corrs[unique_pert] = {}
            for istep,step in enumerate(template_names):
                template1 = step
                if step == 'linear':
                    template_curve = np.array(pert_doses)
                elif step == 'log10':
                    template_curve = np.log10(pert_doses)
                elif step == 'log2':
                    template_curve = np.log2(pert_doses)
                else:
                    print 'template name error'
                # calcualte stats on observation of interest 
                # cc_list = [stats.pearsonr(pert_data[x,:],template_curve) for x in range(len(rids))]
                # rho_vec = [cc_list[x][0] for x in range(len(rids))]
                # rho_vec = np.array(rho_vec)
                # p_vec = [cc_list[x][1] for x in range(len(rids))]
                # p_vec = np.array(p_vec)
                # pvecDictParametric[unique_pert][template1] = p_vec
                # probe_template_corrs[unique_pert][template1] = rho_vec
                # run permutations to creat null distribution of corr values
                ### full matrix of permutations
                nMtrxPerm = n_permutation/len(rids) + 1 #number of matrix permutations needed to reach desired probe perms
                permRhoMtrx = np.zeros((len(rids),nMtrxPerm))
                for perm in range(nMtrxPerm):
                    iRandObs = range(pert_data.shape[1])
                    np.random.shuffle(iRandObs)
                    corrs = np.corrcoef(template_curve,pert_data[:,iRandObs])
                    permRhoMtrx[:,perm] = corrs[0,1:]
                    #test to see if two calculations methods are the same to a given precision
                    # cc_list = [stats.pearsonr(pert_data[x,iRandObs],template_curve) for x in range(len(rids))] #this takes too long
                    # rho_vec = [cc_list[x][0] for x in range(len(rids))]                 
                    # np.allclose(perm_list, np.array(rho_vec),rtol=1e-06)
                #calculate p-value based on null distribution
                grtrMtrx = np.greater(np.abs(permRhoMtrx.T),np.abs(rho_vec))
                null_pVec1 = (1 + np.sum(grtrMtrx, axis=0)) / float(nPerm)
                #compare observed gene to all null genes
                rho_vec = dp.probe_template_corrs[unique_pert][step]
                null_pVec = np.zeros_like(rho_vec)
                for igene in range(len(rho_vec)):
                    rho = rho_vec[igene]
                    p = np.sum(np.abs(permRhoMtrx.flatten()) > np.abs(rho)) /float(len(permRhoMtrx.flatten()))
                    null_pVec[igene] = p
                ## select probe perms1
                num_probes = self.gct.matrix.shape[0]
                probe_inds = range(num_probes)
                perm_cc = []
                for i in range(n_permutation):
                    perm_curve_inds = [random.sample(probe_inds,1)[0] for x in range(len(pert_doses))]
                    perm_curve = [pert_data[perm_curve_inds[x],x] for x in range(len(pert_doses))]
                    perm_covar = np.corrcoef(perm_curve,template_curve)
                    perm_cc.append(perm_covar[0][1])
                corr_null_distribution[unique_pert][template1] = perm_cc
                null_pVec = np.zeros_like(rho_vec)
                for igene in range(len(rho_vec)):
                    rho = rho_vec[igene]
                    p = np.sum(np.abs(perm_cc) > np.abs(rho)) /float(len(perm_cc))
                    null_pVec[igene] = p
                pvecDictEmpirical[unique_pert][template1] = null_pVec
                ### thresholding
                q = .1 #FDR threshold
                pID, pN = FDR.FDR(p_vec,q) #find FDR threshold
                if type(pID) == list:
                    print unique_pert + 'matching to ' + template1 + ' template - perterbation does not have any significant genes that pass the FDR threshold'
                    templateMatchInd[unique_pert][template1] = []
                    continue
                else:
                    pass_fdr = np.less_equal(p_vec,pID) 
                    ipass_fdr = np.array(range(len(rids)))[pass_fdr] #get indices which pass fdr
                    iRhoSort = np.argsort(rho_vec[ipass_fdr])[::-1]
                    iRhoSorted_passFDR = ipass_fdr[iRhoSort] #these are indices which pass FDR and are sorted by correlation
                    data_pass_fdr = pert_data[iRhoSorted_passFDR,:]
                    ordered_rids = [rids[i] for i in iRhoSorted_passFDR]
                    templateMatchInd[unique_pert][template1] = iRhoSorted_passFDR
        self.templateMatchInd = templateMatchInd
        self.pvecDictParametric = pvecDictParametric
        self.probe_template_corrs = probe_template_corrs
        # self.pvecDictEmpirical = pvecDictEmpirical
        # self.corr_null_distribution = corr_null_distribution