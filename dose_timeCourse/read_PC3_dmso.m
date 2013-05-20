fAff = '/xchip/cogs/data/build/affogato/affogato_r1_score_n398050x22268.gctx';
fdmsoLst = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3/6H/wteslm/scratch/dmsoPC3_list.grp';
dmsoLst = importdata(fdmsoLst);

%import LM rids
fASG = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
dbLM = parse_gctx(fASG);
LMrids = dbLM.rid;

%load in LM scores for DMSOs in PC3 cells
db = parse_gctx(fAff,'cid', dmsoLst,'rid',LMrids);

% grab random cid from group DMSO list
r = a + (b-a).*rand(100,1);