#! /usr/bin/env python
'''
create a class/ object to handel common info for dose analysis
'''

import os
import cmap.io.gct as gct
import numpy as np

class DosePlate:
    '''
    class to serve as a containing for dose analysis
    '''
    def __init__(self,s_thresh=6,c_thresh=0.2):
        '''
        constructor
        '''
        self.src = []
        self.c = []
        self.s = []
        self.pid = []
        self.pert_descs = []
        self.pert_ids = []
        self.doses = []
        self.dosePertIDs = []

    def add_from_gct(self,src,ss_column_name='distil_ss', cc_column_name='distil_cc_q75'):
        '''
        reads the meta data of the given gct or gctx file 
        '''
        #set the src for the SC object
        self.src = src
        
        #read in the gct data
        gct_obj = gct.GCT(src=src)
        gct_obj.read()
        
        #grab the pid, ss, and cc data as well as ss and cc cutoffs
        s = gct_obj.get_column_meta(ss_column_name)
        c = gct_obj.get_column_meta(cc_column_name)
        pert_descs = gct_obj.get_column_meta('pert_desc')
        pert_ids = gct_obj.get_column_meta('pert_id')
        doses = gct_obj.get_column_meta('pert_dose')
        id_list = gct_obj.get_column_meta('id')
        pert_desc_list = gct_obj.get_column_meta('pert_desc')
        pid = [x + '::' + pert_desc_list[i] for i,x in enumerate(id_list)]  
        
        #ensure that s and c are lists not numpy arrays
        self.c = list(self.c)
        self.s = list(self.s)
        
        #convert ss and cc into float values
        s = [float(x) for x in s]
        c = [float(x) for x in c]
        
        #add pid, ss, and cc to the existing data
        self.pid.extend(pid)
        self.s.extend(s)
        self.c.extend(c)
        self.pert_ids.extend(pert_ids)
        self.pert_descs.extend(pert_descs)
        self.doses.extend(doses)

    def examine_doses_tested(self):
        '''
        returns:
        doseIndDict - the gct column indices for every unique perturbation on the plate - sorted by dose
        doseSortDict - sorted doses tested for every unique perturbation on the plate
        perts_at_dose - which compounds were tested at more than one dose
        uniqueDoses - unique doses tested if dose testing was consistant across perturbations
        
        *identifies discrepancies in the doses tested for each compound
        '''
        self.uniqueDoses = []
        self.perts_at_dose = []

        qPert = self.pert_descs
        qPertID = self.pert_ids
        qDose = self.doses
        pertSet = set(qPertID)
        #create dictionary of indices, ordered by dose for each compound
        doseIndDict = {}
        doseSortDict = {}
        for pert in pertSet:
            cid_inds = [i for i,x in enumerate(qPertID) if pert in x]
            pert_doses = [float(qDose[x]) for x in cid_inds]
            tmp_tup = zip(pert_doses,cid_inds)
            tmp_tup.sort()
            pert_doses,cid_inds = zip(*tmp_tup)
            doseIndDict[pert] = cid_inds
            doseSortDict[pert] = pert_doses
        self.doseIndDict = doseIndDict
        self.doseSortDict = doseSortDict
        #check to see if compounds were tested with the same number of doses
        compare1 = []
        tested_at_dose = []
        uniqueDoses_list = [] #nested list of doses
        for k,v in self.doseSortDict.iteritems():
            uniqueDoses = set(v)
            nUniqueDoses = len(uniqueDoses)
            if nUniqueDoses < 2:
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
                        print k + ' has a different dose setups on this plate - ' + str(lenComp2) + ' doses vs. ' + str(lenComp1)
        equal_compare = uniqueDoses_list and all(uniqueDoses_list[0] == x for x in uniqueDoses_list) #unique doses are the same for all perts
        if equal_compare == True:
            uniqueDoses = [float(x) for x in uniqueDoses_list[0]]
            uniqueDoses.sort()
            self.uniqueDoses = uniqueDoses
            self.perts_at_dose = tested_at_dose
        else:
            print 'not all perturbations were tested with the same dose set'



    # def clear_points(self):
    #     '''
    #     adds the given ids and indices to the QuadInfo object
    #     '''
    #     self.ids = []
    #     self.indices = []
    #     self.id_map = {}