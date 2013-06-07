%% create a rank file for each CGS gctx

work_dir = '/xchip/cogs/projects/DOS/6June2013b'
Platedir = '/xchip/cogs/projects/DOS/data/pc'
list_dir  =dir([Platedir '/DOSBIO00*'])
for i = 1:length(list_dir);
    PlateName = list_dir(i).name;
    path1 = fullfile(Platedir,PlateName,'by_pert_id_pert_dose');
    path2  =dir([path1 '/' PlateName '*SCORE_LM*']);
    fileName = fullfile(path1,path2(1).name);
    r=regexp(PlateName,'_','split');
    cell = r(2);
    cellList(i) = cell;
    %make gmt file
    outdir = fullfile(work_dir,cell);
    mkdir(outdir{1});
    ds = parse_gctx(fileName)
    [up50, dn50] = get_genesets(ds, 50, 'descend');
    outUP =  fullfile(outdir,[PlateName '_up50.gmt']);
    outDN =  fullfile(outdir,[PlateName '_dn50.gmt']);
%      fid = fopen(outUP{1},'w')
    mkgmt(outUP{1}, up50);
    mkgmt(outDN{1}, dn50);
end

