import cmap.analytics.sc as sc
import cmap.io.gct as gct
import cmap.analytics.dose as doseClass
import cmap.util.progress as progress

# execfile('/xchip/cogs/hogstrom/scripts/dev_pestle/dose_class_obj.py')
#brew by rna_well
gctfile = '/xchip/obelix/pod/brew_tmp/pc/PRISM001_PC3_6H/by_rna_well/PRISM001_PC3_6H_COMPZ.MODZ_SCORE_LM_n260x978.gctx'

dp = doseClass.DosePlate()
dp.add_from_gct(gctfile)
dp.examine_doses_tested()
dp.match_template()
dp.permutation_template()
dp.find_modulated(2,4)
#test ASG file
# ASGfile = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
# ASGdose = doseClass.DosePlate()
# ASGdose.add_from_gct(ASGfile)
# ASGdose.examine_doses_tested()


# sco = sc.SC()
# sco.add_sc_from_gctx_meta(gctfile, verbose=False)
# dose = [float(x.split('::')[0].split(':')[2]) for x in sco.pid]



# p-value by permutation 


'''
matches data to a dose template 
returns p-values, FDR correction 

must add_from_gct() before running
'''
prog = progress.DeterminateProgressBar('template matching')
gcto = dp.gct

doses = gcto.get_column_meta('pert_dose')
perts = gcto.get_column_meta('pert_id')
unique_perts = list(set(perts))
rids = gcto.get_rids()

examineList = dp.perts_at_dose
num_perts = len(examineList)
templateMatchInd = {} #nested dict for index of significant probes
for icmpd,unique_pert in enumerate(examineList):
    prog.update('template match {0}'.format(unique_pert),icmpd,num_perts)
    cid_inds = [i for i,x in enumerate(perts) if unique_pert in x]
    pert_doses = [float(doses[x]) for x in cid_inds]
    tmp_tup = zip(pert_doses,cid_inds)
    tmp_tup.sort()
    pert_doses,cid_inds = zip(*tmp_tup)
    pert_data = gcto.matrix[:,cid_inds]
    if len(pert_doses) > 1:
        template_names = ['linear', 'log10', 'log2']
        templateMatchInd[unique_pert] = {}
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
            cc_list = [stats.pearsonr(pert_data[x,:],template_curve) for x in range(len(rids))]
            rho_vec = [cc_list[x][0] for x in range(len(rids))]
            rho_vec = np.array(rho_vec)
            p_vec = [cc_list[x][1] for x in range(len(rids))]
            p_vec = np.array(p_vec)
            ##run permutations to create a null distribution of corr values
            nPerm = 1000
            permRhoMtrx = np.zeros((len(rids),nPerm))
            for perm in range(nPerm):
                iRandObs = range(pert_data.shape[1])
                np.random.shuffle(iRandObs)
                corrs = np.corrcoef(template_curve,pert_data[:,iRandObs])
                permRhoMtrx[:,perm] = corrs[0,1:]
                #test to see if two calculations methods are the same to a given precision
                # cc_list = [stats.pearsonr(pert_data[x,iRandObs],template_curve) for x in range(len(rids))] #this takes too long
                # rho_vec = [cc_list[x][0] for x in range(len(rids))]					
                # np.allclose(perm_list, np.array(rho_vec),rtol=1e-06)
            ### calculate p-value based on null distribution
            #compare observed gene to null gene 
            grtrMtrx = np.greater(np.abs(permRhoMtrx.T),np.abs(rho_vec))
            null_pVec1 = (1 + np.sum(grtrMtrx, axis=0)) / float(nPerm)
            #compare observed gene to all null genes
            null_pVec = np.zeros_like(rho_vec)
            for igene in range(len(rho_vec)):
                rho = rho_vec[igene]
                p = np.sum(np.abs(permRhoMtrx.flatten()) > np.abs(rho)) /float(len(permRhoMtrx.flatten()))
                null_pVec[igene] = p
    			# compute the random correlation profile for this pert
			
			num_probes = dp.gct.matrix.shape[0]
			probe_inds = range(num_probes)
			perm_cc = []
			# log_perm_cc = []
			# half_log_perm_cc = []
			# quarter_log_perm_cc = []
			for i in range(1000):
				perm_curve_inds = [random.sample(probe_inds,1)[0] for x in range(len(pert_doses))]
				perm_curve = [pert_data[perm_curve_inds[x],x] for x in range(len(pert_doses))]
				# perm_curve = [pert_data[perm_curve_inds[x],x] for x in range(len(pert_doses))]
				perm_covar = np.corrcoef(perm_curve,template_curve)
				perm_cc.append(perm_covar[0][1])
				# log_perm_cc.append(perm_covar[0][2])
				# half_log_perm_cc.append(perm_covar[0][3])
				# quarter_log_perm_cc.append(perm_covar[0][4])
			# compute the nominal p values for all correlation values
			# probe_corrs_p = np.array([stats.percentileofscore(linear_perm_cc,x) 
			# 						for x in rho_vec])
			# p = np.sum(np.abs(perm_cc) > np.abs(rho)) /float(len(perm_cc))
            null_pVec2 = np.zeros_like(rho_vec)
            for igene in range(len(rho_vec)):
                rho = rho_vec[igene]
                p = np.sum(np.abs(perm_cc) > np.abs(rho)) /float(len(perm_cc))
                null_pVec2[igene] = p

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

### edit dose object for new HOG plates
gct1 = '/xchip/obelix/pod/roast/HOG001_MCF7_24H_X1_B10_DUO52HI53LO/zs/HOG001_MCF7_24H_X1_B10_DUO52HI53LO_ZSPCQNORM_n374x978.gct'

dp = doseClass.DosePlate()
dp.add_from_gct(gct1)
dp.examine_doses_tested()
verbose=True

compare1 = []
tested_at_dose = []
uniqueDoses_list = [] #nested list of doses
for k,v in dp.doseSortDict.iteritems():
    uniqueDoses = set(v)
    nUniqueDoses = len(uniqueDoses)
    if nUniqueDoses < 2:
        if verbose:
            print k + ' was tested at only one dose on the plate'
        continue
    else: #save information for perturbations tested at more than one dose
        tested_at_dose.append(k)
        uniqueDoses_list.append(uniqueDoses)
        if not compare1:
            compare1 = uniqueDoses
            lenComp1 = len(compare1)
            continue
        else:
            compare2 = uniqueDoses
            lenComp2 = len(compare2)
            if not compare1 == compare2:
                if verbose:
                    print k + ' has a different dose setups on this plate - ' + str(lenComp2) + ' doses vs. ' + str(lenComp1)
equal_compare = uniqueDoses_list and all(uniqueDoses_list[0] == x for x in uniqueDoses_list) #unique doses are the same for all perts

for x in dp.perts_at_dose:
    print x
    for y in dp.doseSortDict[x]:
        print y
