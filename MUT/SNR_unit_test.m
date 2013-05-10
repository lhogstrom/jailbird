%% load gct filev

    fname = '/xchip/cogs/projects/unit_tests/Diabetes_collapsed_symbols.gct';
    exStruct = parse_gct(fname);
    
    cls_f = '/xchip/cogs/projects/unit_tests/Diabetes.cls';
    clsStruct = importdata(cls_f,' ',2);
    
    rnk_f = '/xchip/cogs/projects/unit_tests/Diabetes.rnk';
    rnkStruct = importdata(rnk_f);
%% Diabetes - run marker selection function 
    
    cls_file = '/xchip/cogs/projects/unit_tests/Diabetes.cls';
    gct_file = '/xchip/cogs/projects/unit_tests/Diabetes_collapsed_symbols.gct';
    outdir = '/xchip/cogs/projects/unit_tests/output_tests/Diabetes';
%     snr_list = markerSelec_SNR(gct_file,outdir,'cls',cls_file,'flip_contrast',true);
    snr_list = markerSelec_SNR(gct_file,outdir,'cls',cls_file);
%     snr_list = markerSelec_SNR(gct_file,outdir,'cls',cls_file,'edge_size',50,'std_fixflow',true);    
    rnk_f = '/xchip/cogs/projects/unit_tests/Diabetes.rnk';
    rnkStruct = importdata(rnk_f);
    exStruct = parse_gct(gct_file);
%% all_aml - run marker selection function 
    
    cls_file = '/xchip/cogs/projects/unit_tests/all_aml_train.cls';
    gct_file = '/xchip/cogs/projects/unit_tests/all_aml_train.gct';
    outdir = '/xchip/cogs/projects/unit_tests/output_tests';
    snr_list = markerSelec_SNR(gct_file,outdir,'cls',cls_file);
    rnk_f = '/xchip/cogs/projects/unit_tests/';
    rnkStruct = importdata(rnk_f);
    exStruct = parse_gct(gct_file);
    
%% p53
    cls_file = '/xchip/cogs/projects/unit_tests/P53.cls';
    gct_file = '/xchip/cogs/projects/unit_tests/P53_collapsed_symbols.gct';
    outdir = '/xchip/cogs/projects/unit_tests/output_tests/p53';
    snr_list = markerSelec_SNR(gct_file,outdir,'cls',cls_file,'std_fixflow',true);
    
    %% load into temp memorry
            exStruct = parse_gct(gct_file);
            clsStruct = importdata(cls_file,' ',2);
            rnk_f = '/xchip/cogs/projects/unit_tests/P53_collapsed_symbols.rnk';
            rnkStruct = importdata(rnk_f);
            
%% run s2n crosscheck 
    %gct_file = '/xchip/cogs/projects/unit_tests/Diabetes_collapsed_symbols.gct';
    %cls_file = '/xchip/cogs/projects/unit_tests/Diabetes.cls';
    cls_file = '/xchip/cogs/projects/unit_tests/P53.cls';
    gct_file = '/xchip/cogs/projects/unit_tests/P53_collapsed_symbols.gct';
    exStruct = parse_gct(gct_file);
    ge = exStruct.mat;
    [CL,CN,NL] = parse_cls(cls_file);
    nclass = 2;
    
    sn = S2N(ge,NL,nclass);

%% make sweet looking 3d heat map

    topEx = qnorm1.mat;
    topEx_mean = mean(topEx')';
    %calculate SD across the row
    topEx_std = std(topEx,0,2);
%     ztopEx = (topEx - topEx_mean) ./ topEx_std;
%     flage = 1;  
    ztopEx = zscore(topEx,1,2);
    bar3(ztopEx)
%         colormap(1)
    
    %whole zscore