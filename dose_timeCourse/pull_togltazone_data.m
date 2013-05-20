fname1 = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_6H/by_pert_id_pert_dose/ASG001_MCF7_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
fname2 = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_24H/by_pert_id_pert_dose/ASG001_MCF7_24H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
fname3 = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
fname4 = '/xchip/obelix/pod/brew/pc/ASG001_PC3_24H/by_pert_id_pert_dose/ASG001_PC3_24H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'

fin = fname1
fout = '/xchip/cogs/projects/ASG_dose_time/data/troglitazone_dose_MCF7_6H.gctx'
db = parse_gctx(fname1)
pertDescs = db.cdesc(:,db.cdict('pert_desc'))
i1 = strmatch('troglitazone',pertDescs)
sigs = db.cid(i1)
db2 = parse_gctx(fname1,'cid',sigs)
mkgctx(fout,db2)
