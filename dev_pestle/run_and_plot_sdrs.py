import cmap.io.gct as gct
import cmap.analytics.dose as doseClass
import cmap.util.progress as progress
import os 

# execfile('/xchip/cogs/hogstrom/scripts/dev_pestle/dose_class_obj.py')
#brew by rna_well
gctfile = '/xchip/obelix/pod/brew_tmp/pc/PRISM001_PC3_6H/by_rna_well/PRISM001_PC3_6H_COMPZ.MODZ_SCORE_LM_n260x978.gctx'

dp = doseClass.DosePlate()
dp.add_from_gct(gctfile)
dp.examine_doses_tested()

cmpd = dp.perts_at_dose[1]
cids = dp.doseIndDict[cmpd] 
rids = dp.gct.get_rids()
AllDoses = dp.gct.get_column_meta('pert_dose')
cmpdDoses = [AllDoses[x] for x in cids]


### build text file for sdrs analysis
work_dir = '/xchip/cogs/hogstrom/analysis/sdrs_perl'
scF = os.path.join(work_dir, 'CMAP-AZD-1152HQPA_dose_data.txt')
headers = list(cmpdDoses) #headers asigned to the dose setup of the last compound
headers.insert(0,'Doses')
with open(scF,'w') as f:
	f.write('\t'.join(headers) + '\n')
	for i,rid in enumerate(rids):
		f.write(rid + '\t')
		for cid in cids: #loop through dose index for the compound
			f.write(str(dp.gct.matrix[i,cid]) + '\t')
		f.write('\n')