#! /usr/bin/env python

'''
-Analyze EMT signatures from Wai Leong Tam et. al. 2013
-run CMAP query on EMT signatures

LH 07/2014
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import cmap.io.gmt as gmt

# set path
wkdir = '/xchip/cogs/hogstrom/bathe/emt_signature'
# geo ID plus differential expression prefix:
signature_dict = {'GSE23952_Sartor':'TGFb_vs_cntrl',
    'GSE17708_Keshamouni':'16h_and_24h_TGBFb_vs_cntrl',
    'GSE43495_Tam': 'TWIST_vs_ctrl_vec',
    'GSE52327_Birnbaum': 'ALDH_plus_vs_ALDH_minus',
    'GSE16838_Vintonenko';'INV_subpopulation_vs_REF',
    'GSE38310_Zhang':'erlotinib_resistant3_DMSO_vs_cntrl_DMSO',
    'GSE38121_Zhang':'erlotonib_resistant_1_2_vs_cntrl'}

### gene symbols --> probes (one-to-many)
### does this artificially inflate genes with more probe ids?
### should it be converted to one-to-one? By choosing one probe randomly? 
def gene_symb_to_probe_id(symbol_list):
    '''
    -given a list of gene symbols, return affy U133A probe ids

    '''
    mc = mu.MongoContainer()
    # geneInfo = mc.gene_info.find({'pr_id':{'$in':list(probe_matrix.index.values)}},
    #         {'pr_id':True,'pr_gene_symbol':True},toDataFrame=True)
    geneInfo = mc.gene_info.find({'pr_gene_symbol':{'$in':symbol_list}},
            {'pr_id':True,'pr_gene_symbol':True},toDataFrame=True)
    geneInfo = geneInfo[~geneInfo.pr_id.isnull()]
    probe_list = list(geneInfo.pr_id.values)
    return probe_list

gmtListUp = []
gmtListDn = []
for sig in signature_dict:
    sig_path = os.path.join(wkdir,sig)
    sigPrefix = signature_dict[sig]
    eFile = os.path.join(sig_path,sigPrefix+'.txt')
    ds = pd.read_csv(eFile,sep='\t')
    ds = ds[~ds['Gene.symbol'].isnull()] # remove rows w/out gene symbols
    ds = ds.sort('t',ascending=False) # sort by t value (positive is up-regulated)
    ds.index = ds['Gene.symbol']
    ds['t_rank'] = ds.t.rank(ascending=False)
    # write up symbols
    outUp = os.path.join(sig_path,sigPrefix+'_100_up_regulated_genes.txt')
    upGenes = ds.ix[:100,'Gene.symbol']
    upGenes.to_csv(outUp,index=False)
    # write dn symbols
    outDn = os.path.join(sig_path,sigPrefix+'_100_dn_regulated_genes.txt')
    dnGenes = ds.ix[-100:,'Gene.symbol']
    dnGenes.to_csv(outDn,index=False)
    # convert to probe ids 
    upProbes = gene_symb_to_probe_id(list(upGenes.values))
    dnProbes = gene_symb_to_probe_id(list(dnGenes.values))
    # save in gmt list 
    # UP
    gmtDictUp = {}
    gmtDictUp['id'] = sig
    gmtDictUp['desc'] = sigPrefix
    gmtDictUp['sig'] = upProbes
    gmtListUp.append(gmtDictUp)
    # Dn
    gmtDictDn = {}
    gmtDictDn['id'] = sig
    gmtDictDn['desc'] = sigPrefix
    gmtDictDn['sig'] = dnProbes
    gmtListDn.append(gmtDictDn)
    # make query directory
    queryDir = os.path.join(sig_path,'cmap_query')
    if not os.path.exists(queryDir):
        os.mkdir(queryDir)
    # write gmt file
    gmtOutUp = os.path.join(queryDir,sig + '_up.gmt')
    gmt.write(gmtListUp,gmtOutUp)
    gmtOutDn = os.path.join(queryDir,sig + '_dn.gmt')
    gmt.write(gmtListDn,gmtOutDn)
    ### run cmap query
    metric = 'wtcs'
    cmd = ' '.join(['rum -q local -f sig_query_tool',
             '--uptag ' + gmtOutUp,
             '--dntag ' + gmtOutDn,
             '--metric ' + metric,
             '--row_space full',
             '--column_space gold',
             '--out ' + queryDir,
             '--mkdir false',
             '--save_tail false'])
    os.system(cmd)
    ### run summly
    # cmd = ' '.join(['rum -q local -f sig_summly_tool',
    #          queryDir,
    #          '--group_query false',
    #          '--out ' + queryDir])
    # os.system(cmd)



# run through MSIG_DB

### EMT questions
# 1) enrichment for any apirori genes
# 2) stability across multiple signatures or TFs
# 3) any drug class enrichment? 

# run up/dn signatures through MSIGDB - do they capture other emt signatures? 


# load kegg pathways
file_kegg = '/xchip/cogs/hogstrom/bathe/gordonov/c2.cp.kegg.v4.0.symbols.gmt'
gt = gmt.read(file_kegg)
keggFrm = pd.DataFrame(gt)
GeneList = keggFrm[keggFrm.id == 'KEGG_REGULATION_OF_ACTIN_CYTOSKELETON'].sig.values
GeneList = list(GeneList[0])

###
aprioriList = ['RAC1',
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
GeneList.extend(aprioriList)
GeneList = list(set(GeneList))


for sig in signature_list:
    eFile = os.path.join(sig_path,sig+'_vs_empty_vec.txt')
    ds = pd.read_csv(eFile,sep='\t')
    ds = ds[~ds['Gene.symbol'].isnull()] # remove rows w/out gene symbols
    ds = ds.sort('t',ascending=False) # sort by t value (positive is up-regulated)
    ds.index = ds['Gene.symbol']
    ds['t_rank'] = ds.t.rank(ascending=False)
    dsMatch = ds[ds['Gene.symbol'].isin(GeneList)]
    x = dsMatch.t_rank 
    plt.hist(x,30,range=(x.min(), x.max()))
    plt.ylabel('gene set members',fontweight='bold')
    plt.xlabel('rank order',fontweight='bold')
    plt.title(sig + ' - actomyosin cytoskeleton gene counts',fontweight='bold')
    outF = os.path.join(sig_path,sig+'_actomyosin_cytoskeleton_enrichment.png')
    plt.savefig(outF, bbox_inches='tight',dpi=200)
    plt.close() 

