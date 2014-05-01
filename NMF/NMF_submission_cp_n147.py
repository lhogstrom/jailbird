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

nComponents = 20
dimDict = {'A375_c20_lm_epsilon':'n1057x978',
'A549_c20_lm_epsilon':'n1405x978', # 
'HA1E_c20_lm_epsilon':'n1345x978',
'HCC515_c20_lm_epsilon':'n1205x978',
'HEPG2_c20_lm_epsilon':'n783x978',
'HT29_c20_lm_epsilon':'n899x978',
'MCF7_c20_lm_epsilon':'n1472x978',
'PC3_c20_lm_epsilon':'n1327x978',
'VCAP_c20_lm_epsilon':'n1359x978'}

#########################
### Define functions ##
#########################

def probe_id_to_gene_symb(inFile,outFile):
    '''
    -change the first column of probe_ids in a gct to gene symbols

    '''
    mtrx = pd.read_csv(inFile,sep='\t',skiprows=[0,1],index_col=0) #,header=True
    probe_ids = mtrx.index.values
    mc = mu.MongoContainer()
    geneInfo = mc.gene_info.find({'pr_id':{'$in':list(probe_ids)}},{'pr_id':True,'pr_gene_symbol':True},toDataFrame=True)
    geneInfo.index = geneInfo.pr_id
    geneInfo = geneInfo.reindex(mtrx.index.values)
    mtrx.index = geneInfo.pr_gene_symbol.values
    mtrx.index.name = 'Name'
    mtrx.to_csv(outFile,sep='\t')
    line_pre_adder(outFile,str(mtrx.shape[0])+'\t'+str(mtrx.shape[1]-1))
    line_pre_adder(outFile,"#1.2")

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

#########################
### Run NMF projection ##
#########################

#specifications for subprocess
processes = set()
max_processes = 9 
### run jobs
wkdir = '/xchip/cogs/projects/NMF/cpd_groups_n147'
for prefix in dimDict:
    print prefix
    dim = dimDict[prefix]
    arg1 = wkdir + '/' + prefix # working directory
    arg2 = 'clique_compound_classes_' + dim
    cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/NMF_code_no_viz.v1.R', # 
         arg1,
         arg2])
    # os.system(cmd)
    processes.add(subprocess.Popen(cmd,shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes if p.poll() is not None)

################################
### Run component annotation ###
################################

#specifications for subprocess
# processes = set()
# max_processes = 9 
# ### run jobs
# outprefix = 'component_annotation_W_rnaseq_projection'
# wkdir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation2'
# for prefix in dimDict:
#     print prefix
#     dim = dimDict[prefix]
#     inDir = wkdir + '/' + prefix # input directory
#     arg1 = inDir + '/' + outprefix # output directory
#     arg2 = 'clique_compound_classes_' + dim
#     arg3 = str(nComponents)
#     ### create output directory
#     if not os.path.exists(arg1):
#         os.mkdir(arg1)
#     ### write gene_symbols to gct
#     inFile = inDir + "/" + arg2 + ".W.k" + str(nComponents) + ".gct"
#     outFile = arg1 + "/" + arg2 + ".W.k" + str(nComponents) + ".gct"
#     probe_id_to_gene_symb(inFile,outFile)
#     cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/Annotate_components_clique.R', # 
#          arg1,
#          arg2,
#          arg3])
#     # os.system(cmd)
#     processes.add(subprocess.Popen(cmd,shell=True))
#     if len(processes) >= max_processes:
#         os.wait()
#         processes.difference_update(
#             p for p in processes if p.poll() is not None)

