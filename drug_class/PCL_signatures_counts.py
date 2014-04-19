'''
-Go through the compounds that make up various PCLs
-Count how many sigatures there are of each compound and the distribution
of cell lines
-identify 'weak' compounds with less signatures 

Larson Hogstrom, 4/2014
'''
import numpy as np
import matplotlib.pyplot as plt
import cmap.util.mongo_utils as mu
import cmap.io.gct as gct
import pandas as pd
import os
import cmap.io.gmt as gmt
import cmap.util.progress as update

# PCL annotation file
# gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140402/cliques.gmt' # gene groupings
# gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140221/cliques.gmt' #most up-to date drug groups
# gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140213/cliques.gmt'
# gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140212/cliques.gmt'
# gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/pcl_20140211/cliques.gmt'
gFile = '/xchip/cogs/projects/pharm_class/rnwork/cliques/cpd_groups_n147.gmt'

wkdir = '/xchip/cogs/projects/pharm_class/lhwork/cpd_groups_n147.gmt_stats'
if not os.path.exists(wkdir):
    os.mkdir(wkdir)

cliqueGMT = gmt.read(gFile)
cliqFrm = pd.DataFrame(cliqueGMT)
# unstack nested list
cliqMemberLong = [item for sublist in cliqFrm.sig.values for item in sublist]
cliqMemb = list(set(cliqMemberLong))

# load summly matrix
summMtrx = '/xchip/cogs/projects/connectivity/summly/matrices/matched_mrp4_n7147x7147.gctx'
gt = gct.GCT()
gt.read(summMtrx)
summFrm = gt.frame

### get info on drug signatures
MC = mu.CMapMongo()
pertInfo = MC.find({'pert_id':{'$in':cliqMemb}},{'sig_id':True,'cell_id':True,'pert_id':True,'pert_iname':True,'is_gold':True},toDataFrame=True)

# tabulate signature stats
pertGrped = pertInfo.groupby('pert_id')
nDrugs = len(pertGrped.groups)
Zs = np.zeros((nDrugs,7))
sig_counts = pd.DataFrame(Zs,index=pertGrped.groups.keys(),columns=['pert_iname','n_cells','n_signatures','n_is_gold','fraction_gold','n_PCLs','PCLs'])
sig_counts.index.name = 'pert_id'
for xt in pertGrped:
    brd = xt[0]
    grp = xt[1]
    sig_counts.ix[brd,'pert_iname'] = grp.pert_iname.values[0]
    cell_set = set(grp.cell_id)
    sig_counts.ix[brd,'n_cells'] = len(cell_set)
    sig_counts.ix[brd,'n_signatures'] = grp.shape[0]
    sig_counts.ix[brd,'n_is_gold'] = sum(grp.is_gold)
    sig_counts.ix[brd,'fraction_gold'] = sum(grp.is_gold)/float(grp.shape[0])
    # check which PCLS the compound belongs to
    def is_in(x):
        return brd in x
    inPCL = cliqFrm[cliqFrm.sig.apply(is_in)]
    sig_counts.ix[brd,'n_PCLs'] = inPCL.shape[0]
    sig_counts.ix[brd,'PCLs'] = list(inPCL.desc.values)
sig_counts['is_summly'] = sig_counts.index.isin(summFrm.index)
outF = os.path.join(wkdir, 'PCL_signature_counts.txt')
sig_counts.to_csv(outF, sep='\t')

plt.hist(sig_counts.n_signatures,30,range=(0,200))
plt.ylabel('number of compounds',fontweight='bold')
plt.xlabel('n_signatures ',fontweight='bold')
plt.title('signature counts for PCL members')
outF = os.path.join(wkdir, 'PCL_n_sigatures.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

plt.hist(sig_counts.n_cells,30)
plt.ylabel('number of compounds',fontweight='bold')
plt.xlabel('n_cells ',fontweight='bold')
plt.title('cell line counts for PCL members')
outF = os.path.join(wkdir, 'PCL_n_cell_lines.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

plt.hist(sig_counts.fraction_gold,30)
plt.ylabel('number of compounds',fontweight='bold')
plt.xlabel('is_gold fraction',fontweight='bold')
plt.title('is_gold signature fraction for PCL members')
outF = os.path.join(wkdir, 'PCL_is_gold_fraction.png')
plt.savefig(outF, bbox_inches='tight',dpi=200)
plt.close()

### write stats on non-summly signatures
non_summly = sig_counts[~sig_counts.is_summly]
outF = os.path.join(wkdir, 'non_summly_PCL_signature_counts.txt')
non_summly.to_csv(outF, sep='\t')


### scratch
# order by group size
cliqFrm['group_size'] = cliqFrm.sig.apply(len)
cliqFrm = cliqFrm.sort('group_size')

for x in cliqFrm['group_size']:
    print x 

gFile = '/xchip/cogs/sig_tools/sig_cliqueselect_tool/sample/PCL_expansion_Apr2014/CID-Synonym-filtered'
cid = pd.read_csv(gFile,sep='\t',header=None,names=['CID','synonym'])
cid[cid.CID == 33]

synList = ['Alexidine dihydrochloride',
'BVT 948',
'NSC 87877',
'Sodium orthovanadate',
'TCS 401',
'Nutlin-3',
'NSC 207895',
'Nutlin-3a',
'Nutlin-3b',
'Caylin-1',
'HLI 373',
'Caylin-2',
'JNJ 26854165',
'NSC 66811',
'MDM2 Inhibitor',
'AS 1892802',
'Fasudil hydrochloride',
'GSK 269962',
'GSK 429286',
'H 1152 dihydrochloride',
'Glycyl-H 1152 dihydrochloride',
'HA 1100 hydrochloride',
'SB 772077B dihydrochloride',
'SR 3677 dihydrochloride',
'Y-27632 dihydrochloride',
'STR1720',
'EX 527',
'Resveratol',
'Sirtinol',
'SIRT1/2 Inhibitor VII',
'Bentamapimod',
'AEG 3482',
'Curcumin',
'BI 78D3',
'CC-401',
'Piceatannol',
'JNK Inhibitor VIII',
'JNK Inhibitor V',
'SP600125',
'SU 3327',
'4-Hydroxynonenal',
'Aloisine A',
'RWJ 67657',
'PKR Inhibitor',
'NU 7441',
'NU 7026',
'PI 103 hydrochloride',
'Compound 401',
'DMNB',
'KU 0060648',
'Doxercalciferol',
'EB 1089',
'Ercalcitriol']



matched = cid[cid.synonym.isin(synList)]
matched.sort('synonym')
