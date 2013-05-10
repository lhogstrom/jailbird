function snr_val = markerSelec_SNR(gct,out,varargin)
% markerSelec_SNR preforme marker selection using signal-to-noise ratio (or
% fold change if there are only one sample for experiment/control)
%   mat_list = markerSelec_SNR(gct,cls_file)
%
%   Inputs:
%           gct
%           out: output directory for ranked list, filtered gctx and heat map
%           class_field1  = column annotation specifying experiment and control 
%           class1_name1 = value of first condition
%           class1_name2 = value of second condition
%           cls - a cls file can be used to specify the SNR contrast
%           gctx file - if gct file is used as an input, this fill will be
%           copied as a gctx file in the input directory
%           flip_contrast - set as 'true' to switch SNR contrast
%           std_fixflow - when set as 'true', the minimum std will be set to 
%           20% of the mean expression value when calculating SNR
%           
%   Outputs:
%           snr_list: snr values for the list edges
%           gene list output: text file containing n genes at the top and
%           bottom of the SNR list
%           heat map: image of heat map
%
% LJH - CMAP, Broad Institute

%% optional parameters
    pnames = {'class_field1','class1_name1','class1_name2','cls','edge_size','flip_contrast','std_fixflow'};
    %edge length
    dflts = {{}, {}, {}, {}, 100,false,false};

    args = parse_args(pnames, dflts, varargin{:});

%% load data
    fname = gct;
    if (strfind(fname,'gctx') >=1) %if gctx load
        mutStruct = parse_gctx(fname);
    else
        sprintf('gct file converting to gctx') %if gct load and write gctx
        mutStruct = parse_gct(fname);
        fnamex = [fname, 'x'];
        mkgctx(fnamex,mutStruct);
        [nr,nc] = size(mutStruct.mat);
        [p, f] = fileparts(fname);
        %strip old dimension if it exists
        prefix = rm_filedim(f);
        gct = fullfile(p, sprintf('%s_n%dx%d.gctx', prefix, nc, nr));
    end
%% class selections
    if ~isempty(args.cls)
        [CL,CN,NL] = parse_cls(args.cls);
        iSect_class1 = find(NL == 1);
        iSect_class2 = find(NL == 2);
    end
    
    if ~isempty(args.class1_name1)
        class_field1 = args.class_field1; 
        class1_name1 = args.class1_name1;
        class1_name2 = args.class1_name2;
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
    end
    
%% flip contrast if applicable
    if args.flip_contrast
        expSet_indx = iSect_class1;
        cntSet_indx = iSect_class2;
    else %use default
        expSet_indx = iSect_class2;
        cntSet_indx = iSect_class1;
    end

%% calculate mean expression for gene 
    bothSet_indx = cat(1,expSet_indx,cntSet_indx);
    
    expSet_mat = mutStruct.mat(:,expSet_indx);
    cntSet_mat = mutStruct.mat(:,cntSet_indx);  
    bothSet_mat = mutStruct.mat(:,bothSet_indx);
    
    %calculate mean across the row if there are more than one rows
        %for exp matrix
        [lnCol, lnRowE] = size(expSet_mat);
        if lnRowE == 1
            expSet_mean = expSet_mat;
        else
            expSet_mean = mean(expSet_mat')';
        end

        %for control matrix
        [lnCol, lnRowC] = size(cntSet_mat);
        if lnRowC == 1
            cntSet_mean = cntSet_mat;
        else
            cntSet_mean = mean(cntSet_mat')';
        end
        
    %% SNR calc or fold changev
    %calculate SD across the row 
    if lnRowE >= 2 & lnRowC >= 2 %skip calculating the standard deviation if there is only one sample
        expSet_std = std(expSet_mat,0,2); 
        cntSet_std = std(cntSet_mat,0,2);    
        expSet_std_fix = std_fixlow(expSet_std, expSet_mean); %make minimum std 20% of the mean
        cntSet_std_fix = std_fixlow(cntSet_std, cntSet_mean); %make minimum std 20% of the mean
       
        % SNR calc or fold changev
                %   (SNR) = (meanExp - meanControl) / (stdExp + stdControl)     
        if args.std_fixflow
            snr_val = (expSet_mean - cntSet_mean) ./ (expSet_std_fix + cntSet_std_fix); %based on std with fixed min    
        else
            snr_val = (expSet_mean - cntSet_mean) ./ (expSet_std + cntSet_std); %based on actual std   
        end
              % replace NaNs with zero
              nan_locations = find(isnan(snr_val));
              snr_val(nan_locations) = 0;
    else
        snr_val = (expSet_mean - cntSet_mean) ./ expSet_mean;
    end
%% SNR calc or fold change

%     %if there is only one sample for experiment/ control - calculate fold
%     %change. Otherwise, calculate SNR
%     if lnRowE ==1 & lnRowC == 1
%         snr_val = (expSet_mean - cntSet_mean) / expSet_mean;
%     else %if there are multiple samples for exp or control, calculate SNR        
%         %   (SNR) = (meanExp - meanControl) / (stdExp + stdControl)     
%         if args.std_fixflow
%             snr_val = (expSet_mean - cntSet_mean) ./ (expSet_std_fix + cntSet_std_fix); %based on std with fixed min    
%         else
%             snr_val = (expSet_mean - cntSet_mean) ./ (expSet_std + cntSet_std); %based on actual std   
%         end
%               % replace NaNs with zero
%               nan_locations = find(isnan(snr_val));
%               snr_val(nan_locations) = 0;
%     end

%% create edge list
    [sort_snr,sort_snr_indx] = sort(snr_val,'descend');
    gene_len = length(mutStruct.mat(:,1)); %number of genes
    edge_size = args.edge_size;
    edge_top = [1:edge_size];
    edge_bottom = [(gene_len-(edge_size-1)):gene_len];
    edges_indx = cat(2,edge_top,edge_bottom)'; 
    
    edges_indx_sort = sort_snr_indx(edges_indx); %the index of genes at the top/bottom  of the snr list
    
    %use gene symbol
%     gene_list = mutStruct.rdesc(edges_indx_sort, mutStruct.rdict('pr_gene_symbol')); % gene symb for top/bottom of list
%     gene_list_top = mutStruct.rdesc(edges_indx_sort(edge_top), mutStruct.rdict('pr_gene_symbol'));
%     gene_list_bottom = mutStruct.rdesc(edges_indx_sort(edge_size+1:(edge_size*2)), mutStruct.rdict('pr_gene_symbol'));
    %use gene rid
    gene_list = mutStruct.rid(edges_indx_sort, 1); % gene symb for top/bottom of list
    gene_list_top = mutStruct.rid(edges_indx_sort(edge_top), 1);
    gene_list_bottom = mutStruct.rid(edges_indx_sort(edge_size+1:(edge_size*2)), 1);
    
    snr_list = snr_val(edges_indx_sort); %snr values for the top and bottom
    mat_list = bothSet_mat(edges_indx_sort,:);
    
%% output files/images
  
    % use Josh's heatmap tool - must be run from the server or have Josh's
    % tools mounted locally 
    
    rid_edges = mutStruct.rid(edges_indx_sort); %define edges to be used
    cid_samples = mutStruct.cid(bothSet_indx);
    rid_grp = [out,'/gene_edges.grp'];
    cid_grp = [out,'/exp_cntl_samples.grp'];
    mkgrp(rid_grp,rid_edges); %make file for filtering - define genes
    mkgrp(cid_grp,cid_samples); %make file for filtering - define samples
    qnorm1 = parse_gctx(gct,'rid',rid_grp,'cid',cid_grp); %make filtered structure 
    filt_name = [out,'/filtered'];
    mkgctx(filt_name,qnorm1); %write filtered structure
    m = length(cid_samples);
    n = length(rid_edges);
    filt_name2 = [out,'/filtered_n',num2str(m),'x',num2str(n),'.gctx']; %set new file name
    
    outfile = [out,'/heatmap_img']; %path and prefix for heatmap image
    [ofname, status, result] = mkheatmap(filt_name2, outfile); %wrapper for josh's heatmap tool
    
    %make grp file - top gene list
    grp_top = [out,'/gene_list_top.grp'];
    mkgrp(grp_top,gene_list_top);
    
    %make grp file - bottom gene list
    grp_bottom = [out,'/gene_list_bottom.grp'];
    mkgrp(grp_bottom,gene_list_bottom);    
    
%     %print top gene list to a file
%         strings = gene_list_top; 
%         fid = fopen(fullfile(out,'gene_list_top.txt'),'w');
%         for row = 1:size(strings,1)
%             fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
%             fprintf(fid, '%s\n', strings{row,end});
%         end
%         fclose(fid);
%         
%     %print bottom gene list to a file 
%         strings = gene_list_bottom; 
%         fid = fopen(fullfile(out,'gene_list_bottom.txt'),'w');
%         for row = 1:size(strings,1)
%             fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
%             fprintf(fid, '%s\n', strings{row,end});
%         end
%         fclose(fid);    
        