%% load gctx

    %fname = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_PC3_6H_X4_B7_DUO52HI53LO_QNORM_n374x978.gct';
    fname1 = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_MCF7_24H_X1_B7_DUO52HI53LO_ZSVCQNORM_n369x978.gct';
    ASG1 = parse_gct(fname1);
    
    fname2 = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_PC3_24H_X2_B7_DUO52HI53LO_ZSVCQNORM_n378x978.gct';
    ASG2 = parse_gct(fname2);   
    
    fname3 = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_PC3_6H_X3_B7_DUO52HI53LO_ZSVCQNORM_n365x978.gct';
    ASG3 = parse_gct(fname3);      
    %plates seem to differ by cell type, time point, replicate
%% view annotations

    dList1 = ASG1.cdesc(:,ASG1.cdict('pert_dose'));
    dListMat1 = cell2mat(dList1); 
    pertList1 = ASG1.cdesc(:,ASG1.cdict('pert_desc'));
    
    pertList2 = ASG2.cdesc(:,ASG2.cdict('pert_desc'));
    pertList3 = ASG3.cdesc(:,ASG3.cdict('pert_desc'));
    
    [uniqCmpds,~,iCmpds] = unique(pertList1);
    
    [C1,IA1,IC1] = unique(pertList1); % compounds tested on plate
    [C2,IA2,IC2] = unique(pertList2);
    [C3,IA3,IC3] = unique(pertList3);
    
%% loop through each unique compound on a plate - identify check number of doses
     
    ASG = ASG1; %plate being examined
    for i = 2:length(uniqCmpds) %skip DSMO (1)
        cmd_indx = find(iCmpds == i); %find indices of each compound
        doseList = ASG.cdesc(cmd_indx,ASG.cdict('pert_dose')); %extract doses for each compound
        doseList = cell2mat(doseList); %convert strings to float
        [unDoseList, ~, idose] = unique(doseList);
        ndoseList(i) = length(unDoseList); %number of unique doses
        
        %determine moderate zscore across all instances of a compound
        zs = ASG.mat(:,cmd_indx); %all the z scores for a given compound
        ridx = 1:size(zs,1); %must give index of zs rows to include
        lmModz = modzs(zs,ridx); %calculate moderated zscore
        
        %sort differental expression for each compound
        [sLM, iLmSort] = sort(lmModz,'descend');
        [~,rank] = sort(iLmSort); %get rank of each gene
        n_edge = 40;
        n_genes = length(lmModz);
        iTop = iLmSort(1:n_edge); %row index of top genes
        iBottom = iLmSort((n_genes-(n_edge-1)):n_genes); %row index of bottom genes
        
        gene1 = iTop(1); %gene of interest
        iTopGeneList(i) = iTop(1);
        n_replicates = 4; %number of replicates for each dose
        zRepMtx = nan(length(unDoseList),n_replicates);
        %loop through each dose
        for j = 1:length(unDoseList)
            iRep = find(idose == j); %find the doses replicates
            iReplicate = cmd_indx(iRep); %index of replicates
            
            %find the zscores of each replicate - one gene
            zRep = ASG.mat(gene1,iReplicate);
            lR = length(zRep); %number of replicates
            zRepMtx(j,1:lR) = zRep; %column = dose, replicates =rows
            
            %for edges, create matrix 
            zRep = ASG.mat(iTop,iReplicate);
            zRepColl = mean(zRep,2); %collapse by averaging zScore
            zCollMtx(j,:) = zRepColl; %dose x gene x replicate
        end
        zClustMtx(i,:,:) = zCollMtx; %compound x dose x gene
        zCmpdMtx(i,:,:) = zRepMtx;
    end
    
    %% plot dose response curve
    c = 18;
    figure(1)
    mtx(:,:) = zCmpdMtx(c,:,:);
        boxplot(mtx',unDoseList,'plotstyle','compact')
        xlabel('dose')
        ylabel('zscore')
        cmpd1 = uniqCmpds(c);
        igene = iTopGeneList(c);
        geneName = ASG.rdesc(igene,ASG.rdict('pr_gene_symbol'));
        titleString = sprintf('compound = %s, gene = %s',cmpd1{1},geneName{1});
        title(titleString)
       
    %graph of dose response for n top genes
    figure(2)
    mtx2(:,:) = zClustMtx(c,:,:);
        boxplot(mtx2',unDoseList,'plotstyle','compact')
        xlabel('dose')
        ylabel('zscore')
        cmpd1 = uniqCmpds(c);
        igene = iTopGeneList(c);
        geneName = ASG.rdesc(igene,ASG.rdict('pr_gene_symbol'));
        titleString2 = sprintf('compound = %s, gene = top 10',cmpd1{1});
        title(titleString2)
    
    
    
     %    
%% sig info
%     sig_info();