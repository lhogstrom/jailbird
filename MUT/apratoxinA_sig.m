%% load in apratoxinA sig
clear
%     fname1 = '/xchip/cogs/hogstrom/analysis/MUT/apratoxinA/aparatoxinA_3h_up.grp';
%     fname2 = '/xchip/cogs/hogstrom/analysis/MUT/apratoxinA/aparatoxinA_3h_dn.grp';
    fname1 = '/xchip/cogs/hogstrom/analysis/MUT/apratoxinA/aparatoxinA_6h_up.grp';
    fname2 = '/xchip/cogs/hogstrom/analysis/MUT/apratoxinA/aparatoxinA_6h_dn.grp';
    aprA_up_list = importdata(fname1);
    aprA_dn_list = importdata(fname2); 
    
%% load chip file

    fid = fopen('/xchip/cogs/hogstrom/analysis/MUT/apratoxinA/chip_id_symbol.txt');
    %fname_chip = '/xchip/cogs/data/vdb/chip/affx.chip';
    % read column headers
    C_text = textscan(fid, '%s %s %s', 'delimiter', '\t');
    pr_id = C_text(1);
    pr_id = pr_id{:};
    geneSym = C_text(3);
    geneSym = geneSym{:};
    
%% match gene symbol in signature file and database
     sig_list = aprA_dn_list;
     [C,ia,ic] = unique(sig_list,'first'); %find repeats 
     start_ndx = (1:size(sig_list))'; %index of starting list
     ndx_repeats = setdiff(start_ndx,sort(ia)); %list repeats
     repeats_list = sig_list(ndx_repeats);
     
     ndx_noRepeate = sort(ia); 
     sig_list = sig_list(ndx_noRepeate); %take out repeats
     ind_lst = [];
     
    %loop through each gene in the signature file 
    for i = 1:length(sig_list)
        apSym = sig_list(i);
        mtch = strmatch(apSym,geneSym,'exact');
        
        if isempty(mtch) %skip if no match
           empty_ndx1(i) = i; 
           
        else
            %eval-(sprintf('Str.%s = mtch;', aprA_dn_list{i}));
            if exist('ind_lst','var')
                ind_lst = cat(1,ind_lst, mtch);
            else
                ind_lst = mtch;
            end
        end
    end
    
    probe_ids = pr_id(ind_lst);
    [g1 empty_ndx2] = find(empty_ndx1 >=1);
    empty_list = sig_list(empty_ndx2); % list out gene Symb with no match
    
        %print gene list to a file
        out = '/xchip/cogs/hogstrom/analysis/MUT/apratoxinA';
        strings = probe_ids; 
        fid = fopen(fullfile(out,'aparatoxinA_6h_dn_prID.grp'),'w');
        for row = 1:size(strings,1)
            fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
            fprintf(fid, '%s\n', strings{row,end});
        end
        fclose(fid);
    
        %print skipped genes to a file
        out = '/xchip/cogs/hogstrom/analysis/MUT/apratoxinA';
        strings = empty_list; 
        fid = fopen(fullfile(out,'aparatoxinA_6h_dn_skipped.txt'),'w');
        for row = 1:size(strings,1)
            fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
            fprintf(fid, '%s\n', strings{row,end});
        end
        fclose(fid);
        
        
        