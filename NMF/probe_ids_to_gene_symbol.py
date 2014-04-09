'''
convert probe ids to gene symbos
'''
import cmap.util.mongo_utils as mu
import pandas as pd

# wFile = '/xchip/cogs/projects/NMF/MCF7_comp_annot_to_CCLE_space/MCF7_top_intra_connecting_compound_classes_n129x978.W.k9.gctx'
# gt = gct.GCT()
# gt.read(wFile)

pFile = '/xchip/cogs/projects/NMF/MCF7_comp_annot_to_CCLE_space/probe_ids.txt'
probes = pd.read_csv(pFile)
probeSer = probes['probe_ids']

mc = mu.MongoContainer()
geneInfo = mc.gene_info.find({'pr_id':{'$in':list(probeSer.values)}},{'pr_id':True,'pr_gene_symbol':True},toDataFrame=True)
geneInfo.index = geneInfo.pr_id
oFile = '/xchip/cogs/projects/NMF/MCF7_comp_annot_to_CCLE_space/probe_id_tbl.txt'
geneInfo = geneInfo.reindex(probeSer.values)
geneInfo.to_csv(oFile,sep="\t")

