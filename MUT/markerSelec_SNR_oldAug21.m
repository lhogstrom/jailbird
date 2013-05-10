function snr_list = markerSelec_SNR(gct_file,class_field1,class1_name1,class1_name2,outdir)
% markerSelec_SNR preforme marker selection using signal-to-noise ratio
%   mat_list = markerSelec_SNR(gct_file,cls_file)
%
%   Inputs:
%           gct_file
%           class_field1  = column annotation specifying experiment and control 
%           class1_name1 = value of first condition
%           class1_name2 = value of second condition
%           outdir: output directory for ranked list and heat map
%   Outputs:
%           snr_list: snr values for the list edges
%           gene list output: text file containing n genes at the top and
%           bottom of the SNR list
%           heat map: image of heat map
%
% LJH - CMAP, Broad Institute
% 8/20/2012

%% load data

    mutStruct = parse_gctx(gct_file);

%% class selections

    indx_annot = mutStruct.cdict(class_field1); %find index in header
    list_annot = mutStruct.cdesc(:,indx_annot); %list annotations for a given header
    
    indx_class1 = strmatch(class1_name1,list_annot); %index of controls
    indx_class2 = strmatch(class1_name2,list_annot); %index of experiment
    
    %match based on a second field if specified (requires matchfield1 and
    %match_name)
    if (exist('match_field1','var'))
        %class matches - make class1 and class2 match for this annotations
        indx_annot2 = mutStruct.cdict(match_field1); %find index in header
        list_annot2 = mutStruct.cdesc(:,indx_annot2); %list annotations for a given header
        indx_matchf1 = strmatch(match_name,list_annot2); %index of experiment

        iSect_class1 = intersect(indx_class1,indx_matchf1);
        iSect_class2 = intersect(indx_class2,indx_matchf1);
    else
        iSect_class1 = indx_class1;
        iSect_class2 = indx_class2;
    end
%% calculate SNR for each gene
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
    %   (SNR) = (meanExp - meanControl) / (stdExp + stdControl)

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
    gene_list_top = mutStruct.rdesc(edges_indx_sort(edge_top), mutStruct.rdict('pr_gene_symbol'));
    gene_list_bottom = mutStruct.rdesc(edges_indx_sort(edge_size+1:(edge_size*2)), mutStruct.rdict('pr_gene_symbol'));
    snr_list = snr_val(edges_indx_sort); %snr values for the top and bottom
    mat_list = bothSet_mat(edges_indx_sort,:);
    
%% output files/images
%     img1 = imagesc(snr_list);
%         title('SNR values')
%         colorbar
%     saveas(img1,fullfile(outdir,'snr_edges.jpg'))
%     
%     img2 = imagesc(mat_list);
%         title('expression heat map')
%         colorbar
%         xlabel('samples/replicates')
%         ylabel('rank order')
%     saveas(img2,fullfile(outdir,'ms_heatmap.jpg'))
    
    % use Josh's heatmap tool
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
    
    %print top gene list to a file
        strings = gene_list_top; 
        fid = fopen(fullfile(outdir,'gene_list_top.txt'),'w');
        for row = 1:size(strings,1)
            fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
            fprintf(fid, '%s\n', strings{row,end});
        end
        fclose(fid);
        
    %print bottom gene list to a file
        strings = gene_list_bottom; 
        fid = fopen(fullfile(outdir,'gene_list_bottom.txt'),'w');
        for row = 1:size(strings,1)
            fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
            fprintf(fid, '%s\n', strings{row,end});
        end
        fclose(fid);    
        