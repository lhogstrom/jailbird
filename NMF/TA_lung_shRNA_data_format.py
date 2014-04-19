'''
-run NMF jobs programatically

Larson Hogstrom, 4/2014
'''
import cmap.io.gct as gct
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.util.progress as update
import glob
import shutil
import subprocess
import fileinput
import cmap.util.mongo_utils as mu

def line_pre_adder(filename,line_to_prepend):
    '''
    -add a new line to a text file

    '''
    f = fileinput.input(filename,inplace=1)
    for xline in f:
        if f.isfirstline():
            print line_to_prepend.rstrip('\r\n') + '\n' + xline,
        else:
            print xline,

### KD data
# # wkdir = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/shRNA_pilot/merged'
# wkdir = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/oe_pilot/'
# # load files
# # mFile= wkdir + '/TA.KD009_KD010_A549_96H_mod_sig_id.txt'
# mFile= wkdir + '/TA.OE006_A549_96H_mod_sig_id.txt'
# mtrx = pd.read_csv(mFile,sep='\t',index_col=0)
# aFile = wkdir + '/TA.KD009_KD010_A549_96H_annotations.txt'
# # aFile = wkdir + '/TA.KD009_KD010_A549_96H_cgs_w_plate_info_annotations.txt'
# annt = pd.read_csv(aFile,sep='\t')

### OE data
# wkdir = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/oe_pilot'
# load files
# mFile= wkdir + '/TA.OE006_A549_96H_mod_sig_id.txt'
# mtrx = pd.read_csv(mFile,sep='\t',index_col=0)
# aFile = wkdir + '/TA.OE006_A549_96H_mod_sig_id_annotations.txt'
# annt = pd.read_csv(aFile,sep='\t')

### OE data - w/ replicates
wkdir = '/xchip/cogs/projects/NMF/TA_lung_cancer_genomic_perturbations/oe_pilot/oe_plate_reps'
# load files
mFile= wkdir + '/A549_mini_plate_reps_modHeader.txt'
mtrx = pd.read_csv(mFile,sep='\t',index_col=0)
aFile = wkdir + '/A549_mini_plate_reps_annotations.txt'
annt = pd.read_csv(aFile,sep='\t')



### prep matrix into gct form
# outFile = wkdir + '/TA.KD009_KD010_A549_96H_mod_sig_id.gct'
outFile= wkdir + '/A549_mini_plate_reps_annotations.gct'
mtrx.index.name = 'id'
mtrx.to_csv(outFile,sep='\t')
### add lines for gct headers
line_pre_adder(outFile,str(mtrx.shape[0])+'\t'+str(mtrx.shape[1]-1))
line_pre_adder(outFile,"#1.2")

### make gmts of gene shRNAs
geneGrped = annt.groupby('pert_id')
gmtList = []
for grp in geneGrped:
    gmtDictUp = {}
    gmtDictUp['id'] = grp[0]
    gmtDictUp['desc'] = grp[0]
    gmtDictUp['sig'] = list(grp[1].sig_id.values)
    gmtList.append(gmtDictUp)
# gmtOut = wkdir + '/gene_shRNA_sig_id.gmt'
gmtOut = wkdir + '/gene_oe_sig_id.gmt'
gmt.write(gmtList,gmtOut)

### load core drivers - save sig_ids to new gmt
gFile= wkdir + '/core_lung_drivers.gmt'
coreGMT = gmt.read(gFile)
coreOE = coreGMT['sig']
coreFrm = annt[annt.pert_id.isin(coreOE)]
sig_ids = list(coreFrm.sig_id.values)
gmtDict = {}
gmtDict['id'] = 'core_lung_drivers'
gmtDict['desc'] = 'core_lung_drivers'
gmtDict['sig'] = sig_ids
gmtOut = wkdir + '/core_lung_drivers_sig_id.gmt'
gmt.write([gmtDict],gmtOut)

