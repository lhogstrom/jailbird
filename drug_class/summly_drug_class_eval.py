import os
import cmap.io.gct as gct
import numpy as np
import cmap.util.mongo_utils as mu
import glob
import matplotlib.pyplot as plt
import cmap.analytics.dgo as dgo
import cmap.util.progress as progress
import pandas as pd
import matplotlib
import pandas as pd

work_dir = '/xchip/cogs/hogstrom/analysis/summly/cp_class'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

# drugFile = '/xchip/cogs/projects/target_id/ctd2_annots/ctd2_merged_mapped_genes.txt'
drugFile = '/xchip/cogs/projects/cp_annot/drug_classes_AS.txt'
drugLabels = pd.io.parsers.read_csv(drugFile,sep='\t')
# targetSet = set(drugLabels['gene_dw_annot'])
# labelSet = set(drugLabels['CLASS'])

grpSet = set(drugLabels['CLASS'])
grpToCp = {}
for grp in grpSet:
    grpPerts = drugLabels['BRD'][drugLabels['CLASS'] == grp]
    grpToCp[grp] = list(grpPerts.values)
# compound to group dict
cpToGrp = {}
for ibrd, brd in enumerate(drugLabels['BRD']):
    cpToGrp[brd] = drugLabels['CLASS'][ibrd]

alldrug = drugLabels['BRD']
allinames = drugLabels['PERT_INAME']
pDescDict = {}
for ibrd,brd in enumerate(alldrug):
    pDescDict[brd] = allinames[ibrd]

CM = mu.CMapMongo()
# look up brds for unkown
for ibrd, brd in enumerate(drugLabels['BRD']):
    if brd == '-666':
        iname = drugLabels['PERT_INAME'][ibrd]
        #look up the brd using iname
        trueBRD = CM.find({'pert_iname':iname},{'pert_id':True},limit=1)
        if trueBRD:
            trueBRD = trueBRD[0]
            #replace in dataframe
            drugLabels['BRD'][ibrd] = trueBRD

### look up sig_ids for each brd
drugs = 