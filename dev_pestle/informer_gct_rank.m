%% create a rank file for each CGS gctx




work_dir = '/xchip/cogs/hogstrom/analysis/informer_CTD'
CGSdir = '/xchip/cogs/projects/rnai_analysis/collapsed_signatures/Collapsed_Signatures_Landmark_978'
list_dir  =dir([CGSdir '/*_*_modz_by_pert_desc_signatures_any_target_n*x978.gctx'])
for i = 1:length(list_dir);
    path = fullfile(CGSdir,list_dir(i).name);
    fileName = list_dir(i).name;
    r=regexp(fileName,'_','split');
    cellList(i) = r(1);
    mkdir(work_dir,r{1});
    copyfile(path,fullfile(work_dir,r{1}));
    rnkStruc = score2rank(path);
    rankOut = fullfile(work_dir,r{1},sprintf('%s_modz_CGZ_RANK', r{1}));
    mkgctx(rankOut,rnkStruc)
end

