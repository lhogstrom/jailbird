'''
Find data related to Simon's thesis project

5/15/2014
'''
import pandas as pd
import os
import cmap.io.gmt as gmt
import glob
import shutil
import cmap.util.mongo_utils as mu

### actomyosin cytoskeleton organization
# Keywords:
# Rac
# Cdc42
# Rho-ROCK
# PKC
#  delta-myosin 
# light chain kinase-myosin
# leading edge proteins (Arp 2/3, PI3K)
# MAPK pathways (p38 in particular)
# proteases that direct trailing end retraction such as Calpains
# focal adhesion kinase
# Src

## gene symbols
RAC1
CDC42
RHOA
ROCK1
RICS
RHOA
PRKCA
PIK3CA
ARPC1A
MAPK
ERK
MAPK14
CAPN4
CAPN1
CAPN2
PTK2
SRC
# RHOA pathway
NgR1
LINGO1
p75
TROY
# myosin genes:
MYH3
MYH6
MYH7
MYH9
MYH11
MYO1A
MYO5A
MYO6 
MYO7A
MYO15A

###
GeneList = ['RAC1',
'CDC42', 
'RHOA',
'ROCK1',
'RICS',
'RHOA',
'PRKCA',
'PIK3CA',
'ARPC1A',
'MAPK',
'ERK',
'MAPK14',
'CAPN4',
'CAPN1',
'CAPN2',
'PTK2',
'SRC',
'NgR1',
'LINGO1',
'p75',
'TROY',
'MYH3',
'MYH6',
'MYH7',
'MYH9',
'MYH11',
'MYO1A',
'MYO5A',
'MYO6 ',
'MYO7A',
'MYO15A']

# load kegg pathways
file_kegg = '/xchip/cogs/hogstrom/bathe/gordonov/c2.cp.kegg.v4.0.symbols.gmt'
gt = gmt.read(file_kegg)
keggFrm = pd.DataFrame(gt)
GeneList = keggFrm[keggFrm.id == 'KEGG_REGULATION_OF_ACTIN_CYTOSKELETON'].sig.values
GeneList = list(GeneList[0])

wkdir = '/xchip/cogs/hogstrom/bathe/gordonov'
### genomic perturbation
# shRNA
mc = mu.MongoContainer()
cgsFrm = mc.sig_info.find({'pert_type':'trt_sh.cgs','pert_iname':{'$in':list(GeneList)}},
            {'dn100_full':False,'up100_full':False,'dn100_bing':False,'up100_bing':False,'dn50_lm':False,'up50_lm':False,'pert_idose':False,'pert_dose_unit':False},toDataFrame=True)
cgsFrm.index = cgsFrm.sig_id
outF = wkdir + '/cytoskeleton_shRNA_signatures.txt'
cgsFrm.to_csv(outF, sep='\t')

# OE
mc = mu.MongoContainer()
oeFrm = mc.sig_info.find({'pert_type':'trt_oe','pert_iname':{'$in':GeneList}},
            {'dn100_full':False,'up100_full':False,'dn100_bing':False,'up100_bing':False,'dn50_lm':False,'up50_lm':False,'pert_idose':False,'pert_dose_unit':False},toDataFrame=True)
outF = wkdir + '/cytoskeleton_oe_signatures.txt'
oeFrm.index = oeFrm.sig_id
oeFrm.to_csv(outF, sep='\t')



