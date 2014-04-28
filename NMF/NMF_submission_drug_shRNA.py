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

nComponents = 9
shRNA_dimDict = {'A375_shRNA_c9_lm_epsilon':'n279x978',
'A549_shRNA_c9_lm_epsilon':'n258x978', # 
'HA1E_shRNA_c9_lm_epsilon':'n267x978',
'HCC515_shRNA_c9_lm_epsilon':'n232x978',
'HEPG2_shRNA_c9_lm_epsilon':'n252x978',
'HT29_shRNA_c9_lm_epsilon':'n259x978',
'MCF7_shRNA_c9_lm_epsilon':'n400x978',
'PC3_shRNA_c9_lm_epsilon':'n441x978',
'VCAP_shRNA_c9_lm_epsilon':'n448x978'}

cp_dimDict = {'A375_drug_c9_lm_epsilon':'n340x978',
'A549_drug_c9_lm_epsilon':'n345x978', # 
'HA1E_drug_c9_lm_epsilon':'n380x978',
'HCC515_drug_c9_lm_epsilon':'n364x978',
'HEPG2_drug_c9_lm_epsilon':'n268x978',
'HT29_drug_c9_lm_epsilon':'n301x978',
'MCF7_drug_c9_lm_epsilon':'n450x978',
'PC3_drug_c9_lm_epsilon':'n428x978',
'VCAP_drug_c9_lm_epsilon':'n380x978'}

shRNASer = pd.Series(shRNA_dimDict)
cpSer = pd.Series(cp_dimDict)

cellList = ['A375','A549', 'HA1E', 'HCC515', 'HEPG2', 'HT29', 'MCF7', 'PC3', 'VCAP'] # cmap 'core' cell lines
# make data frames of the meta data
cpFrm = pd.DataFrame(data={'cp_dim':cpSer.values,
                    'cp_name':cpSer.index,
                    'shRNA_dim':shRNASer.values,
                    'shRNA_name':shRNASer.index},
                    index=cellList)

# wkdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA'
# for prefix in dimDict:
#     print prefix
#     dim = dimDict[prefix]
#     inDir = wkdir + '/' + prefix # input directory
#     g = os.listdir(inDir)
#     print g

#move input files from one directory to another
# wkdir = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation'
# wkdir2 = '/xchip/cogs/projects/NMF/NMF_parameter_evaluation2'
# for prefix in dimDict:
#     print prefix
#     dim = dimDict[prefix]
#     mvDir = wkdir2 + '/' + prefix
#     if not os.path.exists(mvDir):
#         os.mkdir(mvDir)
#     gfile = '/' + prefix + '/clique_compound_classes_' + dim + '.gct'
#     aFile = '/' +  prefix + '/clique_compound_classes.v2.txt'
#     shutil.copy(wkdir+gfile, wkdir2+gfile)
#     shutil.copy(wkdir+aFile, wkdir2+aFile)

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

## cp 
# specifications for subprocess
# processes = set()
# max_processes = 9 
# ### run jobs
# wkdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA'
# for prefix in cp_dimDict:
#     print prefix
#     dim = cp_dimDict[prefix]
#     arg1 = wkdir + '/' + prefix # working directory
#     arg2 = 'clique_compound_classes_' + dim
#     cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/NMF_code_no_viz.v1.R', # 
#          arg1,
#          arg2])
#     # os.system(cmd)
#     processes.add(subprocess.Popen(cmd,shell=True))
#     if len(processes) >= max_processes:
#         os.wait()
#         processes.difference_update(
#             p for p in processes if p.poll() is not None)

## shRNA
# # specifications for subprocess
# processes = set()
# max_processes = 9 
# ### run jobs
# wkdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA'
# for prefix in shRNA_dimDict:
#     print prefix
#     dim = shRNA_dimDict[prefix]
#     arg1 = wkdir + '/' + prefix # working directory
#     arg2 = 'shRNA_drug_target_genes_' + dim
#     cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/NMF_code_no_viz.v1.R', # 
#          arg1,
#          arg2])
#     # os.system(cmd)
#     processes.add(subprocess.Popen(cmd,shell=True))
#     if len(processes) >= max_processes:
#         os.wait()
#         processes.difference_update(
#             p for p in processes if p.poll() is not None)

#################################
### drug --> shRNA projections ##
#################################

# processes = set()
# max_processes = 9 
# ### run jobs
# wkdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA'
# for cell in cpFrm.index:
#     # this W matrix will be used
#     path1 = wkdir + '/' + cpFrm.ix[cell,'cp_name']
#     prefix1 = 'clique_compound_classes_' + cpFrm.ix[cell,'cp_dim']
#     # new H matrix will match the above W matrix
#     path2 = wkdir + '/' + cpFrm.ix[cell,'shRNA_name']
#     prefix2 = 'shRNA_drug_target_genes_' + cpFrm.ix[cell,'shRNA_dim']   
#     nFolder = 'drug_projections'
#     nFPath = path2 + '/' +nFolder
#     if not os.path.exists(nFPath):
#         os.mkdir(nFPath)
#     cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/Comp_TA.v3_Drug_shRNA.R', # 
#          path1,
#          prefix1,
#          path2,
#          prefix2,
#          nFolder])
#     # os.system(cmd)
#     processes.add(subprocess.Popen(cmd,shell=True))
#     if len(processes) >= max_processes:
#         os.wait()
#         processes.difference_update(
#             p for p in processes if p.poll() is not None)

#################################
### shRNA --> drug projections ##
#################################

processes = set()
max_processes = 9 
### run jobs
wkdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA'
for cell in cpFrm.index:
    # this W matrix will be used
    path1 = wkdir + '/' + cpFrm.ix[cell,'shRNA_name']
    prefix1 = 'shRNA_drug_target_genes_' + cpFrm.ix[cell,'shRNA_dim']   
    # new H matrix will match the above W matrix
    path2 = wkdir + '/' + cpFrm.ix[cell,'cp_name']
    prefix2 = 'clique_compound_classes_' + cpFrm.ix[cell,'cp_dim']
    nFolder = 'shRNA_projections'
    nFPath = path2 + '/' +nFolder
    if not os.path.exists(nFPath):
        os.mkdir(nFPath)
    cmd = ' '.join(['Rscript /xchip/cogs/hogstrom/scripts/jailbird/NMF/Comp_TA.v3_Drug_shRNA.R', # 
         path1,
         prefix1,
         path2,
         prefix2,
         nFolder])
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
# wkdir = '/xchip/cogs/projects/NMF/NMF_drug_shRNA'
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

