%% load T2D data

    %fname = '/xchip/cogs/stacks/deprecated/STK010_T2D/STK010_QNORM_n756x978.gct';
    fname = '/xchip/cogs/data/tnroast/plates/MUT001_HEK293T_XH_X1_B1_DUO53HI52LO/MUT001_HEK293T_XH_X1_B1_DUO53HI52LO_QNORM_n333x978.gct';
    mutStruct = parse_gct(fname);

%% view contents of annotations

    for i = 1:length(mutStruct.chd)
       annot = mutStruct.chd(i)
       tmp_annot_list = mutStruct.cdesc(:,i);
       
       %convert
       if isnumeric(tmp_annot_list{1}) == 1
           tmp_annot_list = cell2mat(tmp_annot_list);
           unique(tmp_annot_list)
       else
           unique(tmp_annot_list)           
       end
       
    end
         
%% class selections

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
%% find replicates
clear('samelist','difflist')
    for j = 1:length(mutStruct.chd) %loop through each annotation, see if they are the same
        tmp_annot_list = mutStruct.cdesc(indx_class1,j);
        
        %first_copy_list(1:length(tmp_annot_list)) = tmp_annot_list(1); %copy the first item multiple times
        % check to see if all of the array is the same as the first item
           if isnumeric(tmp_annot_list{1}) == 1
               tmp_annot_list = cell2mat(tmp_annot_list);
               first_match = strmatch(tmp_annot_list(1),tmp_annot_list);
           else  
               first_match = strmatch(tmp_annot_list(1),tmp_annot_list);
           end
        
        if length(first_match) == length(tmp_annot_list)
            samelist(j) = j;
        else
            difflist(j) = j;
        end
    end
    
    samelist = samelist(samelist~=0);
    difflist = difflist(difflist~=0);
    
    % see where annotations with the same pert_desc are different
        for i = 1:length(difflist)
           annot = mutStruct.chd(difflist(i))
           value = mutStruct.cdesc(indx_class1,difflist(i))
        end
        
%% find the controls which match up with each experiment
% 
%     qc_annot_indx = [1     2     8     9    12    28]; %fields that will awlays be different (QC or plate locations)
%     %'qc_f_logp'     'qc_iqr'     'count_cv'     'count_mean'     'det_well'     'rna_well' 
%     
%     variable_indx = 1:length(mutStruct.chd);
%     variable_indx(qc_annot_indx) = []; %remove qc annotations from index of interest
%     
%     %compare single control where does/doesn't it match exp
%     k = 10;
%     cntl_match = indx_class1(k);
%     
%     
%     clear('samelist','difflist')
%     for j = 1:length(variable_indx) %loop through each annotation, see if they are the same
%         tmp_annot_list = mutStruct.cdesc(indx_class2,variable_indx(j));
%         
%         cntrl_annot = mutStruct.cdesc{cntl_match,variable_indx(j)}; %pick one cntl experiemnt, find its annotation
%         %first_copy_list(1:length(tmp_annot_list)) = tmp_annot_list(1); %copy the first item multiple times
%         % check to see if all of the array is the same as the first item
%            if isnumeric(tmp_annot_list{1}) == 1
%                tmp_annot_list = cell2mat(tmp_annot_list);
%                first_match = strmatch(cntrl_annot,tmp_annot_list);
%            else  
%                first_match = strmatch(cntrl_annot,tmp_annot_list);
%            end
%         
%         if length(first_match) == length(tmp_annot_list)
%             samelist(j) = variable_indx(j);
%         else
%             difflist(j) = variable_indx(j);
%         end
%     end
    
  %% find the controls witch match up with each experiment

    qc_annot_indx = [1     2     3     8     9    12    28]; %fields that will awlays be different (QC or plate locations)
    %'qc_f_logp'     'qc_iqr'     'count_cv'     'count_mean'     'det_well'     'rna_well' 
    
    variable_indx = 1:length(mutStruct.chd);
    variable_indx(qc_annot_indx) = []; %remove qc annotations from index of interest
    
    %what are the exp conditions
    indx_exp = indx_class2;
    indx_cntl_long = indx_class1; %list of all control samples
   
    clear('samelist','difflist')
    for j = 1:length(variable_indx) %loop through each annotation, see if they are the same
        tmp_annot_list_cntl = mutStruct.cdesc(indx_cntl_long,variable_indx(j));
        tmp_annot_list_exp = mutStruct.cdesc(indx_class2,variable_indx(j));
        
        k = 1; %which of the experimental indices do you want to use as a refernce?
        exp_match = indx_class2(k);
        exp_annot = mutStruct.cdesc{exp_match,variable_indx(j)}; %pick one cntl experiemnt, find its annotation
        %first_copy_list(1:length(tmp_annot_list)) = tmp_annot_list(1); %copy the first item multiple times
        % check to see if all of the array is the same as the first item
           if isnumeric(tmp_annot_list_cntl{1}) == 1
               tmp_annot_list_cntl = cell2mat(tmp_annot_list_cntl);
               first_match = strmatch(exp_annot,tmp_annot_list_cntl);
           else  
               first_match = strmatch(exp_annot,tmp_annot_list_cntl);
           end
        
        if length(first_match) == length(tmp_annot_list_cntl)
            samelist(j) = variable_indx(j);
        else
            difflist(j) = variable_indx(j);
        end
    end  
    
    %check which ones are the same/diff
    annot1 = 15;
    annot_nm = mutStruct.chd(annot1)
    exp_out = mutStruct.cdesc(indx_class2(1),annot1)
    cntl_out = mutStruct.cdesc(indx_cntl_long(1),annot1)


    
%% calculate SNR for each gene
%   (SNR) = (?Exp ? ?Control) / (?Exp+ ?Control)
    %? = mean expression
    %? = SD of expression
    
    expSet_indx = indx_class2;
    cntSet_indx = indx_cntl_long(1:5);
    bothSet_indx = cat(1,expSet_indx,cntSet_indx);
    
    expSet_mat = mutStruct.mat(:,expSet_indx);
    cntSet_mat = mutStruct.mat(:,cntSet_indx);  
    bothSet_mat = mutStruct.mat(:,bothSet_indx);
    
    %calculate mean across the row
    expSet_mean = mean(expSet_mat')';
    cntSet_mean = mean(cntSet_mat')';
   
    %% whats up with this std?
    %calculate SD across the row
        %expSet_std = std(expSet_mat')';
        %cntSet_std = std(cntSet_mat')';    
    expSet_std = std(expSet_mat,0,2);
    cntSet_std = std(cntSet_mat,0,2);    
    %c0sigma = std(c0ge,0,2);
    
    %% 
    %(SNR) = (?Exp ? ?Control) / (?Exp+ ?Control)
    snr_val = (expSet_mean - cntSet_mean) ./ (expSet_std + cntSet_std);
    
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
