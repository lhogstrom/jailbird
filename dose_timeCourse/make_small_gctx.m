%write in a big gctx file, select one type of compound and write the subset
%of data to a gctx file

gctIn = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_24H/by_pert_id_pert_dose/ASG001_MCF7_24H_COMPZ.MODZ_SCORE_LM_n85x978.gctx'
db = parse_gctx(gctIn);
pertDescs = db.cdesc(:,db.cdict('pert_desc'));
matchDrug ='withaferin-a';
keepBin = strcmpi(pertDescs,matchDrug);
ikeep = find(keepBin >= 1);

%add second compound
ikeep(end+1:(end+length(ikeep))) = cellstrfind(pertDescs,'trichostatin-a'); %add to ikeep

cidKeep = db.cid(ikeep);
db2 = parse_gctx(gctIn,'cid',cidKeep);
outdir = '/xchip/cogs/hogstrom/analysis/ASG';
% fname = sprintf('ASG_MCF7_24H_%s.gctx',matchDrug); %use matchDrug to name file
fname = 'ASG_MCF7_24H_smallData.gctx'
dsfile = fullfile(outdir,fname);

mkgctx(dsfile,db2);