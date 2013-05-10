%% load 3 MUT replicate gcts   

%     merge_folder_dataset('plate', {'MUT001_HEK293T_XH_X1_B1_DUO53HI52LO', 'MUT001_HEK293T_XH_X2_B1_DUO53HI52LO', 'MUT001_HEK293T_XH_X3_B1_DUO53HI52LO'}, 'plate_path', '/xchip/cogs/data/tnroast/plates', 'dstype', 'QNORM', 'out', '/xchip/cogs/hogstrom/analysis/MUT', 'use_gctx', false)

%% load MUT data

    %fname = '/xchip/cogs/stacks/deprecated/STK010_T2D/STK010_QNORM_n756x978.gct';
    %fname = '/xchip/cogs/data/tnroast/plates/MUT001_HEK293T_XH_X1_B1_DUO53HI52LO/MUT001_HEK293T_XH_X1_B1_DUO53HI52LO_QNORM_n333x978.gct';
    fname = '/xchip/cogs/hogstrom/analysis/MUT/merged_QNORM_n993x978.gct';
    mutStruct = parse_gct(fname);   
    
    %write gctx file
    fnamex = [fname, 'x'];
    mkgctx(fnamex,mutStruct);
    
%% class selections
    cls_file = '/xchip/cogs/hogstrom/analysis/MUT/lrrk2_class.cls';
    clsStruc = importdata(cls_file,' ');
    
    class_field1 = 'pert_desc'; %header of interest
    indx_annot = mutStruct.cdict(class_field1); %find index in header
    list_annot = mutStruct.cdesc(:,indx_annot); %list annotations for a given header
    
    class1_name = 'LRRK2'; 
    class2_name = 'Empty Vector'; % 'ATG16L1' 'Foxp3' 'LRRK2' 'PTEN'
    indx_class1 = strmatch(class1_name,list_annot); %index of controls
    indx_class2 = strmatch(class2_name,list_annot); %index of experiment
    
    %class matches - make class1 and class2 match for this annotations
    match_field1 = 'pert_2_id';
    indx_annot2 = mutStruct.cdict(match_field1); %find index in header
    list_annot2 = mutStruct.cdesc(:,indx_annot2); %list annotations for a given header
    match_name = 'TB';
    indx_matchf1 = strmatch(match_name,list_annot2); %index of experiment
    
    iSect_class1 = intersect(indx_class1,indx_matchf1);
    iSect_class2 = intersect(indx_class2,indx_matchf1);
    
%% calculate SNR for each gene
%   (SNR) = (?Exp ? ?Control) / (?Exp+ ?Control)
    %? = mean expression
    %? = SD of expression
    
    expSet_indx = iSect_class1;
    cntSet_indx = iSect_class2;
    bothSet_indx = cat(1,expSet_indx,cntSet_indx);
    
    expSet_mat = mutStruct.mat(:,expSet_indx);
    cntSet_mat = mutStruct.mat(:,cntSet_indx);  
    bothSet_mat = mutStruct.mat(:,bothSet_indx);
    
    %calculate mean across the row
    expSet_mean = mean(expSet_mat')';
    cntSet_mean = mean(cntSet_mat')';
   
    %calculate SD across the row
        %expSet_std = std(expSet_mat')';
        %cntSet_std = std(cntSet_mat')';    
    expSet_std = std(expSet_mat,0,2);
    cntSet_std = std(cntSet_mat,0,2);    
    expSet_std_fix = std_fixlow(expSet_std, expSet_mean); %make minimum std 20% of the mean
    cntSet_std_fix = std_fixlow(cntSet_std, cntSet_mean); %make minimum std 20% of the mean
    %% SNR calc
    %snr_val = (expSet_mean - cntSet_mean) ./ (expSet_std + cntSet_std); %based on actual std
    snr_val = (expSet_mean - cntSet_mean) ./ (expSet_std_fix + cntSet_std_fix); %based on std with fixed min
    
    %% create edge list
    [sort_snr,sort_snr_indx] = sort(snr_val,'descend');
    gene_len = length(mutStruct.mat(:,1)); %number of genes
    edge_size = 50;
    edge_top = [1:edge_size];
    edge_bottom = [(gene_len-(edge_size-1)):gene_len];
    edges_indx = cat(2,edge_top,edge_bottom)'; 
    
    edges_indx_sort = sort_snr_indx(edges_indx); %the index of genes at the top/bottom  of the snr list
    
    gene_list = mutStruct.rdesc(edges_indx_sort, mutStruct.rdict('pr_gene_symbol')); % gene symb for top/bottom of list
    snr_list = snr_val(edges_indx_sort); %snr values for the top and bottom
    mat_list = bothSet_mat(edges_indx_sort,:);

%% outputs
    img1 = imagesc(snr_list);
        title('SNR values')
        colorbar
    saveas(img1,'snr_edges.jpg')
    
    img2 = imagesc(mat_list);
        title('expression heat map')
        %colorbar
        t = colorbar('peer',gca); get(t,'ylabel');
        set(t,'ylabel','My Title');
        xlabel('samples/replicates')
        ylabel('rank order')
    saveas(img2,'ms_heatmap.jpg')
    
    strings = gene_list; 
    fid = fopen('gene_list.txt','w');
    for row = 1:size(strings,1)
        fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
        fprintf(fid, '%s\n', strings{row,end});
    end
    fclose(fid);
    %% run marker selection function 
    
%     fname = '/xchip/cogs/data/tnroast/plates/MUT001_HEK293T_XH_X1_B1_DUO53HI52LO/MUT001_HEK293T_XH_X1_B1_DUO53HI52LO_QNORM_n333x978.gct';
    cls_file = '/xchip/cogs/hogstrom/analysis/MUT/lrrk2_class.cls';
    gct_file = '/xchip/cogs/hogstrom/analysis/MUT/merged_QNORM_n993x978.gctx';
    outdir = '/xchip/cogs/hogstrom/analysis/MUT/markerSelection_test';
%     snr_list = markerSelec_SNR(gct_file,cls_file,outdir);
%     snr_list = markerSelec_SNR(gct_file,outdir,'pert_desc','PTEN','Empty Vector');
    snr_list = markerSelec_SNR(gct_file,outdir,'class_field1','pert_desc','class1_name1','PTEN','class1_name2','Empty Vector');
    
%% use Josh's heatmap tool
    rid_edges = mutStruct.rid(edges_indx_sort);
    cid_samples = mutStruct.cid(bothSet_indx);
    rid_grp = [outdir,'/gene_edges.grp'];
    cid_grp = [outdir,'/exp_cntl_samples.grp'];
    mkgrp(rid_grp,rid_edges);
    mkgrp(cid_grp,cid_samples);
    qnorm1 = parse_gctx(gct_file,'rid',rid_grp,'cid',cid_grp);  
    filt_name = [outdir,'/filtered'];
    mkgctx(filt_name,qnorm1);
    m = length(cid_samples);
    n = length(rid_edges);
    filt_name2 = [outdir,'/filtered_n',num2str(m),'x',num2str(n),'gctx'];
    
    outfile = [outdir,'/heatmap_img'];
    [ofname, status, result] = mkheatmap(filt_name2, outfile);