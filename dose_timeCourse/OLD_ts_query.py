cd /xchip/cogs/hogstrom/scripts/dose_timeCourse

	# Given a signature find its attributes
	# pert = signature.find({'sig_id': 'OSG001_SALE_XH:A19'})
	# pert = pert[0]
	# pert.keys()
	
	# 'T2D_L1000analysis_20120912.txt'
	# Short 'T2D_L1000analysis_20120912_short.txt'
	# First Pass 'T2D_L1000analysis_20120912_first_pass.txt'
	#infile = open('T2D_L1000analysis_20120912_first_pass.txt', 'r')
	
	#genes = []
	#for line in infile:
	#line = line.strip("\r")
	#line = line.split("\t")
	#genes.append(line[0])
	#genes.remove('Gene symbol')
	#genes = set(genes)
	
	#infile.close()

import pymongo

conn =  pymongo.Connection('vitalstatistix')
affogato = conn['affogato']
affogato.authenticate('cmap_user', 'l1000')
signature = affogato['signature']

# Make a .grp file with just the signature id, keep only those signatures that pass the following criteria:
# is gold (is_gold > 0),
# signature strength (distil_ss > 4)
#
#
# 'T2Dgenes.grp'
# Short 'T2Dgenes_short.grp'
# First Pass 'T2Dgenes_first_pass.grp'
with open('T2Dgenes_first_pass.grp', 'w') as outfile:
    for gene in genes:
        perts = signature.find({'pert_desc': gene}, {'is_gold': True, 'sig_id': True, 'distil_ss': True})
        res = [pert for pert in perts]
	for r in res:
#	    if r['is_gold'] > 0:
#	        if r['distil_ss'] > 4:
	            outfile.write(str(r['sig_id']) + '\n')
		     
# No Filtering
with open('T2Dgenes_first_pass.grp', 'w') as outfile:
    for gene in genes:
        perts = signature.find({'pert_desc': gene}, {'is_gold': True, 'sig_id': True, 'distil_ss': True})
        res = [pert for pert in perts]
	for r in res:
	    outfile.write(str(r['sig_id']) + '\n')

# Make a .txt file with the gene, cell id, is gold, signature id, and signature strength
# 'T2Dgenes.txt'
# Short 'T2Dgenes_short.txt'
with open('T2Dgenes_first_pass.txt', 'w') as outfile:
     outfile.write('Gene' + '\t' + 'Cell id' + '\t' + 'Is gold' + '\t' + 'Sig id' + '\t' + 'Signature Strength' + '\n')
     for gene in genes:
        perts = signature.find({'pert_desc': gene}, {'cell_id':True,  'is_gold':True, 'sig_id':True, 'distil_ss': True})
        res = [pert for pert in perts]
	if len(res) == 0:
	    outfile.write(gene + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')
	else:
	    for r in res:
	        outfile.write(gene + '\t' + str(r['cell_id']) + '\t' + str(r['is_gold']) + '\t' + str(r['sig_id']) + '\t' + str(r['distil_ss']) + '\n')