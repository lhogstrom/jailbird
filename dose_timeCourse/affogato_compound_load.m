%% list unique compounds for ASG plates
%full
%cmpdList = {'alvespimycin'; 'troglitazone'; 'trichostatin-a'; 'geldanamycin'; 'radicicol'; 'alpha-estradiol'; 'tanespimycin'; 'wortmannin'; 'sirolimus'; 'thioridazine'; 'fulvestrant'; 'estradiol'; 'mitoxantrone'; 'LY-294002'; 'valproic-acid'; 'genistein'; 'fluphenazine'; 'tretinoin'; 'vorinostat'; 'withaferin-a'; 'trifluoperazine'};

%short - take out compounds with instances with mu
cmpdList = {'alvespimycin'; 'troglitazone'; 'trichostatin-a'; 'geldanamycin'; 'radicicol'; 'alpha-estradiol'; 'tanespimycin'; 'wortmannin'; 'thioridazine'; 'fulvestrant'; 'estradiol'; 'valproic-acid'; 'genistein'; 'fluphenazine'; 'tretinoin'; 'vorinostat'; 'withaferin-a'; 'trifluoperazine'};

for i = 1:length(cmpdList)
    cmpd = cmpdList{i};
    f1 = sprintf('%s_byID.grp',cmpd);
    dir1 = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/cmpd_sigIDs';
    fname = fullfile(dir1,f1);
%% load in sig ID list
    %fname = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/cmpd_sigIDs/radicicol_byID.grp';
    sigList = importdata(fname);
    %sigList = sigList(1:1000);
    
%% load in expression data
    if length(sigList) >= 1 %analyze if not empty
        bigGctx = '/xchip/cogs/data/build/affogato/affogato_r1_ranklm_n398050x978.gctx'; %LM set
        db = parse_gctx(bigGctx,'cid',sigList);
        annt = sig_info(sigList); % load annotations for data

        pertList = extractfield(annt,'pert_desc')';
        sigList2 = extractfield(annt,'sig_id')';
        pertIdList = extractfield(annt,'pert_id')';
        doseList = extractfield(annt,'pert_dose')';
        timeList = extractfield(annt,'pert_time')';

        doseTrack(i) = length(unique(doseList));
        pertTrack(i) = length(pertList);
    
    else
        doseTrack(i) = NaN;
        pertTrack(i) = NaN;        
    end
end