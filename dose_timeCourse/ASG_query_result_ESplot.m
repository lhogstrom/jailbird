%% load in relevant compounds and doses    
    %fname = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_n85x22268.gctx';
    %fname = '/xchip/obelix/pod/brew/pc/ASG001_PC3_24H/by_pert_id_pert_dose/ASG001_PC3_24H_COMPZ.MODZ_SCORE_n85x22268.gctx';
    %fname = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_6H/by_pert_id_pert_dose/ASG001_MCF7_6H_COMPZ.MODZ_SCORE_n85x22268.gctx'; %mcf7, 6h
    fname = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_24H/by_pert_id_pert_dose/ASG001_MCF7_24H_COMPZ.MODZ_SCORE_n85x22268.gctx';
    %fname = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_PC3_6H_COMPZ.MODZ_SCORE_n85x22268.gctx';
    %fname = '/xchip/cogs/projects/ASG_dose_time/data/ASG001_MCF7_24H_COMPZ.MODZ_SCORE_n85x22268.gctx';
    db = parse_gctx(fname);
    % set info about setup
    cellLine = 'MCF7';
    timeP = '24H';
        %% 22k space queries
            %specify output of query_tool: 
            %MCF7 - 24H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/nov02/my_analysis.query_tool.6279027';
            %PC3 - 6H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/nov05/my_analysis.query_tool.6727508';
            %PC3 - 24H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/nov05/my_analysis.query_tool.6727542';
        %% lm - all of affogato
            %MCF7 6H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/nov07/my_analysis.query_tool.7297360';
            %MCF7 24H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/nov07/my_analysis.query_tool.7297449';
            %PC3 - 6H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/nov06/my_analysis.query_tool.6953388';
            %PC3 - 24H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/nov08/my_analysis.query_tool.7412638';
        %% lm - only is gold
            %MCF7 6H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/6H_isgld/nov12/my_analysis.query_tool.2012111212201291';
            %MCF7 24H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/24H_isgld/nov12/my_analysis.query_tool.2012111212234791';
            %PC3 - 6H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/6H_isgld/nov09/my_analysis.query_tool.7731750';
            %PC3 - 24H
            %qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/PC3_sigs/LM_queries/24H_isgld/nov09/my_analysis.query_tool.2012110919262691';
        %% results to cross check 
            qdir = '/xchip/cogs/hogstrom/analysis/scratch/Nov20/dose_analysis_tool.1353449771597/nov20/my_analysis.query_tool.2012112017162991';

%% load result gctx
    gctFname3 = fullfile(qdir,'result_ESLM.COMBINED_n85x398050.gctx'); %bigger file for all of affogato
    %gctFname3 = fullfile(qdir,'result_ESLM.COMBINED_n85x198733.gctx'); %smaller is just is_gold instances
    rdb = parse_gctx(gctFname3);
    
    fname1 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_desc__affogato.txt';
        %affPert = importdata(fname1);
        fid =fopen(fname1,'r'); % fid is the file identifier
        affPert=textscan(fid,'%s','delimiter', '');
        fclose(fid)
        affPert=affPert{:};
    fname2 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/pert_id_affogato.txt';
        %affPertID = importdata(fname2);
        fid =fopen(fname2,'r'); % fid is the file identifier
        affPertID=textscan(fid,'%s','delimiter', '');
        fclose(fid)
        affPertID=affPertID{:};
    fname3 = '/xchip/cogs/hogstrom/scripts/dose_timeCourse/sig_ids_affogato.txt';
        % affSigID = importdata(fname3);
        fid =fopen(fname3,'r'); % fid is the file identifier
        affSigID =textscan(fid,'%s','delimiter', '');
        fclose(fid)
        affSigID =affSigID{:};
%% match annotations for result matrix
    rSigID = rdb.rid;
    [sA, iS] = sort(affSigID); %sort affogato sigIDs
    [sR, iR] = sort(rSigID); %sort results sigIDs
    if ~isequal(sA,sR)
        sprintf('conflict between result and affogato sigIDs')
    end
    
    saffPertID = affPertID(iS);
    saffPert = affPert(iS);
    qPertID = db.cdesc(:,db.cdict('pert_id')); %pertIDs for each query
    qPert = db.cdesc(:,db.cdict('pert_desc')); %pertIDs for each query
    
    %resort result matrix to match sigID sort
    ESmat = rdb.mat(iR,:);
    [sESmat, iES] = sort(ESmat,1,'descend'); %sort acording to ES
    
    %outdir for new figure jpgs
    outdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_24H_selfConnect_counts/viz';
    %loop through each of the queries
    for i = 1:length(qPertID)
        IDsorted = saffPertID(iES(:,i));
        qStr = qPertID{i}; %define pertID string
        cmpd1 = db.cdesc(i,db.cdict('pert_desc'));
        dose1 = db.cdesc(i,db.cdict('pert_dose')); 
        if length(qStr) >= 13
            qStr = qStr(1:13); %shorten pertID string
        end
        matchList = strmatch(qStr,IDsorted);
        
%         %hist of ES rank axis 
%             j = figure;
%             hist(matchList,50)
%                 xlabel('ES rank')
%                 ylabel('freq')
%                 outName2 = sprintf('%s_%sUM_%s_%s_EShist.jpg',cmpd1{1},num2str(dose1{1}),cellLine,timeP);
%                 fout = fullfile(outdir,outName2);
%             r = 150; % pixels per inch
%             set(j, 'PaperUnits', 'inches', 'PaperPosition', [0 0 800 200]/r); %dim of image
%             print(j,'-dpng',sprintf('-r%d',r), fout);
%             close(j);
%             
%         %rank plot
%             n_cmpds = length(rSigID);
%             h = figure;
%             hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
%                          'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
%                          'XLim',[1 n_cmpds],...               %#   set the x axis limit,
%                          'YLim',[0 eps],...
%                          'Color','none');               %#   and don't use a background color
%             plot(matchList,0,'b.','MarkerSize',5);  %# Plot data set 2
%             strTitle = rdb.cid{i}
%             title(strrep(strTitle,'_','\_'))
%             outName = sprintf('%s_%sUM_%s_%s_ESscat.jpg',cmpd1{1},num2str(dose1{1}),cellLine,timeP);
%             fout = fullfile(outdir,outName);
%             %control output image
%                 r = 150; % pixels per inch
%                 set(h, 'PaperUnits', 'inches', 'PaperPosition', [0 0 800 200]/r); %dim of image
%                 print(h,'-dpng',sprintf('-r%d',r), fout);
%                 %saveas(h,fout)
%                 hold off;
%             close(h);
    end
    