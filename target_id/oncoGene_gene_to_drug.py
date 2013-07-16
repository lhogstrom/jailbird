import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd

##############################
### make drug name dictionary ###
##############################

### CTD2
work_dir = '/xchip/cogs/projects/target_id/CTD2_25June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

targetSheetF = '/xchip/cogs/projects/target_id/4June2013/Informer2_drug_targets.txt'
targetDict = {}
pDescDict = {}
with open(targetSheetF,'rt') as f:
    for string in f:
        splt = string.split('\r')
        for i,line in enumerate(splt):
            splt2 = line.split('\t')
            pID = splt2[0] #the pert_id listed the line
            pDesc = splt2[1]
            targets = splt2[2]
            targets = targets.split(';')
            targets = [x for x in targets if x != '']
            if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
                continue
            else:
                targetDict[pID] = targets
                pDescDict[pID] = pDesc
### drugbank
work_dir = '/xchip/cogs/projects/target_id/DrugBank_26June2013'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

targetSheetF = '/xchip/cogs/projects/target_id/7June2014/A2_DrugBank_targets_tab.txt'
targetDict = {}
pDescDict = {}
with open(targetSheetF,'rt') as f:
    for string in f:
        splt = string.split('\r')
        for i,line in enumerate(splt):
            splt2 = line.split('\t')
            pID = splt2[0] #the pert_id listed the line
            pDesc = splt2[1]
            targets = splt2[2:]
            targets = [x for x in targets if x != '']
            # targets = targets.split(';')
            if targets[0] == '' or targets[0] == '?' or targets[0] == '-666':
                continue
            else:
                targetDict[pID] = targets
                pDescDict[pID] = pDesc
### DOS
work_dir = '/xchip/cogs/projects/DOS/8July_target_id'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

##############################
### run dgo operations ###
##############################

reload(dgo)
dg = dgo.QueryTargetAnalysis(out=work_dir + '/drug_KD_spearman')
# dg.add_dictionary(targetDict=targetDictCGS)
# dg.get_sig_ids(genomic_pert='KD',is_gold=True)
# dg.run_drug_gene_query(metric='spearman',max_processes=10)
# #wait until queries finish
dg.make_result_frames(gp_type='KD',metric='spearman')

# make empty pDescDict
fullBRDs = []
for ind in dg.dfCS.index:
    brd = ind[0]
    fullBRDs.append(brd)
uniqBRDs = list(set(fullBRDs))
pDescDict = {}
for brd in uniqBRDs:
    pDescDict[brd] = '-666'

inFile = '/xchip/cogs/projects/oncoDome/OncoDome_genes.txt'
outDir = 'oncoDome_11July'
if not os.path.exists(dg.outputdir+'/'+outDir):
    os.mkdir(dg.outputdir+'/'+outDir)
cgsList = []
with open(inFile,'rt') as f:
    for string in f:
        splt = string[:-1]
        cgsList.append(splt)
geneAll = set(cgsList)
# check to see gene has a CGS
geneCGS = geneAll.copy()
CM = mu.CMapMongo()
for gene in geneAll:
	CGSsigs = CM.find({'pert_type':'trt_sh.cgs','pert_iname':gene},{'sig_id':True,'pert_iname':True})
	if not CGSsigs:
		geneCGS.remove(gene)
for gene in geneCGS:
    dg.gene_to_drug_similarity(testGene=gene,
                                gp_type='KD',
                                metric='spearman',
                                outName=outDir + '/gene_to_drug',
                                pDescDict=pDescDict,
                                n_rand=10000,
                                n_uncorrected=20,
                                connection_test='two_sided')

#pax8

##############################
### drugbank xml annotations ###
##############################

### retrieve info from drugbank - gene target info and 
#number of drugs targeting that drug
file1 = '/xchip/cogs/projects/target_id/drugBank_db_xml/drugbank.xml'
context = ET.iterparse(file1, events=("start", "end"))
context = iter(context)
event, root = context.next()
for child in root:
    print child.tag, child.attrib

contextFull = ET.iterparse(file1)
for event, elem in ET.iterparse(file1):
	print event, elem
elem.tag
elem.text


#make list of drugbank drugs
nameList = []
for event, elem in ET.iterparse(file1):
	if elem.tag.endswith("drug"): #"name"
		nameList.append(elem.text)




def parseXML():
    file = open(str(options.drugxml),'r')
    data = file.read()
    file.close()
    dom = parseString(data)
    druglist = dom.getElementsByTagName('drug')

    with codecs.open(str(options.csvdata),'w','utf-8') as csvout, open('DrugTargetRel.csv','w') as dtout:
        for entry in druglist:
        count = count + 1
        try:
            drugtype = entry.attributes['type'].value
            print count
        except:
            print count
            print entry
            drugidObj = entry.getElementsByTagName('drugbank-id')[0]
            drugid = drugidObj.childNodes[0].nodeValue
            drugnameObj = entry.getElementsByTagName('name')[0]
            drugname = drugnameObj.childNodes[0].nodeValue

            targetlist = entry.getElementsByTagName('target')
            for target in targetlist:
                targetid = target.attributes['partner'].value
                dtout.write((','.join((drugid,targetid)))+'\n')

            csvout.write((','.join((drugid,drugname,drugtype)))+'\n')

#!python
import codecs
from xml.dom import minidom

class DrugBank(object):
    def __init__(self, filename):
        self.fp = open(filename, 'r')

    def __iter__(self):
        return self

    def next(self):
        state = 0

        while True:
            line = self.fp.readline()

            if state == 0:
                if line.strip().startswith('<drug '):
                    lines = [line]
                    state = 1
                    continue

                if line.strip() == '</drugs>':
                    self.fp.close()
                    raise StopIteration()

            if state == 1:
                lines.append(line)
                if line.strip() == '</drug>':
                    return minidom.parseString("".join(lines))

with codecs.open('csvout.csv', 'w', 'utf-8') as csvout, open('dtout.csv', 'w') as dtout:
    db = DrugBank('drugbank.xml')
    for dom in db:
	        entry = dom.firstChild
	        drugtype = entry.attributes['type'].value
	        drugidObj = entry.getElementsByTagName('drugbank-id')[0]
	        drugid = drugidObj.childNodes[0].nodeValue
	        drugnameObj = entry.getElementsByTagName('name')[0]
	        drugname = drugnameObj.childNodes[0].nodeValue
	        targetlist = entry.getElementsByTagName('target')
	        for target in targetlist:
	            targetid = target.attributes['partner'].value
            dtout.write((','.join((drugid,targetid)))+'\n')

        csvout.write((','.join((drugid,drugname,drugtype)))+'\n')
