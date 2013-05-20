%% load gctx
    %MCF7 24 h
%     fname1 = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_24H_X1/by_pert_id_pert_dose/ASG001_MCF7_24H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
%     fname2 = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_24H_X4/by_pert_id_pert_dose/ASG001_MCF7_24H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
    fname1 = '/Users/hogstrom/Documents/work_scratch/ASG001_MCF7_24H_COMPZ.MODZ_SCORE_LM_X1_n85x978.gctx';
    fname2 = '/Users/hogstrom/Documents/work_scratch/ASG001_MCF7_24H_COMPZ.MODZ_SCORE_LM_X4_n85x978.gctx';
    
    
    %MCF7 6h
%     fname3 = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_6H_X1/by_pert_id_pert_dose/ASG001_MCF7_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
%     fname4 = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_6H_X2/by_pert_id_pert_dose/ASG001_MCF7_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
    fname3 = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_MCF7_6H_COMPZ.MODZ_SCORE_LM_X1_n85x978.gctx';
    fname4 = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_MCF7_6H_COMPZ.MODZ_SCORE_LM_X2_n85x978.gctx';
    
    
    %PC3 6h
%     fname5 = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H_X3/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
%     fname6 = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H_X4/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
    fname5 = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_X3_n85x978.gctx';
    fname6 = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_X4_n85x978.gctx';    
    
    %import data
    ASG1 = parse_gctx(fname3);
    ASG2 = parse_gctx(fname4);         
    cellLine = 'MCF7';
    timeP = '6H';
    
    %plates seem to differ by cell type, time point, replicate
%% view annotations
    %is check to see if a sepcific field is the same
    isequal(ASG1.cdesc(:,ASG1.cdict('pert_dose')),ASG2.cdesc(:,ASG1.cdict('pert_dose')))

    dList1 = ASG1.cdesc(:,ASG1.cdict('pert_dose'));
    dListMat1 = cell2mat(dList1); 
    pertList1 = ASG1.cdesc(:,ASG1.cdict('pert_desc'));
    pertList2 = ASG2.cdesc(:,ASG2.cdict('pert_desc'));
    
    [uniqCmpds,~,iCmpds] = unique(pertList1);
    
    [C1,IA1,IC1] = unique(pertList1); % compounds tested on plate
    [C2,IA2,IC2] = unique(pertList2);
    %[C3,IA3,IC3] = unique(pertList3);
    
%% loop through each unique compound on a plate - identify check number of dose

    ASG = ASG1; %plate being examined
    SS = sig_strength(ASG.mat); % calculate signature strength 
    doseAll = ASG.cdesc(:,ASG.cdict('pert_dose'));
    
    n_edge = 40; 
    n_doses = 4; %number of dose concentrations on plate
    n_replicates = 1; %number of replicates for each dose
    cmd_indices = nan(length(uniqCmpds),n_doses);
    zClustMtx = nan(length(uniqCmpds),n_doses,n_edge);
    % loop through each unique compound on a plate
    doseAll= cell2mat(doseAll); %convert strings to float 
    for i = 2:length(uniqCmpds) %skip DSMO (1)
        cmd_indxS = find(iCmpds == i); %find indices of each compound
            doseList = ASG.cdesc(cmd_indxS,ASG.cdict('pert_dose')); %extract doses for each compound
            doseList = cell2mat(doseList); %convert strings to float     
            [sd, id] = sort(doseList); %sort the doses for each compound
            cmd_indx = cmd_indxS(id); %re-order cmd_indx so to it goes from lowest to highest dose
        cmd_indices(i,:) = cmd_indx; 

        [unDoseList, ~, idose] = unique(doseList)
        ndoseList(i) = length(unDoseList); %number of unique doses
        
        %determine moderate zscore across all instances of a compound
        zs = ASG.mat(:,cmd_indx); %all the z scores for a given compound
        ridx = 1:size(zs,1); %must give index of zs rows to include
        lmModz = modzs(zs,ridx); %calculate moderated zscore
        
        %sort differental expression for each compound
        [sLM, iLmSort] = sort(lmModz,'descend');
        [~,rank] = sort(iLmSort); %get rank of each gene
        n_genes = length(lmModz);
        iTop = iLmSort(1:n_edge); %row index of top genes
        iBottom = iLmSort((n_genes-(n_edge-1)):n_genes); %row index of bottom genes
        
        gene1 = iTop(1); %gene of interest - over/under expressed gene
        %gene1 = iBottom(n_edge);
        iSingleEdgeGeneList(i) = gene1; %save index edge gene to list
        zRepMtx = nan(length(unDoseList),n_replicates);
        zCollMtx = nan(length(unDoseList),length(iTop));
        %loop through each dose
        for j = 1:length(unDoseList)
            iRep = find(idose == j); %find the doses replicates
            iReplicate = cmd_indxS(iRep); %index of replicates
                        
            %find the zscores of each replicate - one gene
            zRep = ASG.mat(gene1,iReplicate);
            lR = length(zRep); %number of replicates
            zRepMtx(j,1:lR) = zRep; %column = dose, replicates =rows
            
            %for edges, create matrix 
            zRep = ASG.mat(iTop,iReplicate); %look at top genes
            %zRep = ASG.mat(iBottom,iReplicate); %look at top genes
            zRepColl = mean(zRep,2); %collapse by averaging zScore
            zCollMtx(j,:) = zRepColl; %dose x gene x replicate
        end
        zClustMtx(i,:,:) = zCollMtx; %z-score of edges - compound x dose x geneEdge
        %zClustMtx(i,:) = zCollMtx;
        zCmpdMtx(i,:,:) = zRepMtx; %z-score of top gene - compound x dose x replicate
        %zCmpdMtx(i,:) = zRepMtx;
    end
    
%% how does signature strength correlate with dose response?
    figure(3)
    hist(SS,30)
        title('signature strength across all doses')
        xlabel('signature strength')
        ylabel('freq')
    
    figure(4)
    plot(doseAll,SS,'o')    
        title('signature strength by dose')
        xlabel('dose concentraiton')
        ylabel('signature strength')    
    %try some sort of log transformation here?
    
%% plot dose response curve
    %compounds with potential dose response:
    %MCF7 24H- dose response - 20,4,6,9,11,13,15,16 - top
    %MCF7 24H- dose response - 6, 11, 16, 17, 20  - bottom
    
    %MCF7 6H- dose response -  - top
    %MCF7 6H- dose response - 4,5,7,10,11,16,20,21  - bottom    
        
    %PC3 6H- dose response -  - top
    %PC3 6H- dose response - 2, 4(qq plot predicts diff, but not SS), 8, 10, 11, 12,13,16,20,21,22   - bottom
    
    dir1 = sprintf('%s_sigs',cellLine);
    dir2 = sprintf('matlab_plots_%s',timeP);
    outdir = fullfile('/xchip/cogs/projects/ASG_dose_time/cmap_queries/',dir1,dir2);
    mkdir(outdir)
    

    
    c = 20;
    ic = cmd_indices(c,:);
    ASG.cdesc(ic,ASG.cdict('pert_dose')) %print out order of doses
    
    % graph dose response of a single gene (at very edge)
    figure(1)
    mtx(:,:) = zCmpdMtx(c,:,:);
        boxplot(mtx',unDoseList,'plotstyle','compact')
        xlabel('dose')
        ylabel('zscore')
        cmpd1 = uniqCmpds(c);
        igene = iSingleEdgeGeneList(c);
%         igene = iBottom (c);
        geneName = ASG.rdesc(igene,ASG.rdict('pr_gene_symbol'));
        titleString = sprintf('compound = %s, gene = %s',cmpd1{1},geneName{1});
        title(titleString)
       
    %graph of dose response for n top genes
    h = figure;
    mtx2(:,:) = zClustMtx(c,:,:);
        boxplot(mtx2',unDoseList,'plotstyle','compact')
        xlabel('dose')
        ylabel('zscore')
        cmpd1 = uniqCmpds(c);
        geneName = ASG.rdesc(igene,ASG.rdict('pr_gene_symbol'));
        titleString2 = sprintf('compound = %s, genes = top %s',cmpd1{1},num2str(n_edge));
        title(titleString2)
        outName = sprintf('%s_top_genes.jpg',cmpd1{1});
        fout = fullfile(outdir,outName);
        saveas(h,fout)
        
    % qq plot for the experiment
     figure(3)
     zlist1 = sort(ASG.mat(:,ic));
     qqplot(zlist1)    
     hleg1 = legend('d=.08','d=.4','d=2','d=10');
     
     %plot sig strength
     figure(4)
     ss1 = SS(ic);
        boxplot(ss1',unDoseList,'plotstyle','compact','outliersize',10)
        %plot(ss1,'.','MarkerSize',40)
        xlabel('dose')
        ylabel('sig strength')
        title(cmpd1)
        %set(gca,'xticklabel',unDoseList)

%% visualize DMSO qq plot
    % qq plot for the experiment
    iDmso = 1; %dmso is the first column in the data
    dDMSO = ASG1.mat(:,1); %data for the DSMOs
     figure(1)
     qqplot(dDMSO)    
     %hleg1 = legend('d=.08','d=.4','d=2','d=10');

%% qq plot of gene expression
     %histogram of all expression values
     %hist(ASG.mat(:),30) %these tails are very big - does this represent the bioactivity? 
     
    %build normal distribution which reflects the null distribution of z-scores
     nullD = random('Normal',0,1,1,978);
     figure(3)
     zlist1 = sort(ASG.mat(:,ic));
%      qqplot(zlist1(:,1),nullD) 
%      hold on 
%      qqplot(zlist1(:,2),nullD)
%      qqplot(zlist1(:,3),nullD)
     qqplot(zlist1(:,4),nullD)
%      hold off
     
     figure(2)
     qqplot(zlist1,nullD) 
     hleg1 = legend('d=.08','d=.4','d=2','d=10');
     
%      figure(1)
%      qqplot(zlist1) 
%      hleg1 = legend('cos_x','sin_x');
     
%% sig info
%     sig_info();