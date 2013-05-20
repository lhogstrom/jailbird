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
            qdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs/LM_queries/nov07/my_analysis.query_tool.7297449';
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


%% loop through each compound/dose
    pertList = db.cdesc(:,db.cdict('pert_desc')); 
    pertIDList = db.cdesc(:,db.cdict('pert_id')); 
    doseList = db.cdesc(:,db.cdict('pert_dose'));
    for i = 1:length(pertList);
    %for i = 1:4;  
        %i=2
        cmpd = pertList{i};
        cmpdID = pertIDList{i}; 
        if length(cmpdID) >= 13 %check to see if the pertID is less than 13 char (eg DMSO)
                cmpdID_shrt = cmpdID(1:13); %shorten pertID to only 13 char
            else
                cmpdID_shrt = cmpdID; %leave pertID if less than 13 char
        end
        cmpdU = upper(pertList{i});
        dose1 = db.cdesc(i,db.cdict('pert_dose')); 
        fname = sprintf('tail_ESLM_%s_%sUM_%s_%s.txt',cmpdU,num2str(dose1{1}),cellLine,timeP); %name of summary table for each compound
        tblList{i} = fname; 
        
        dirFname = fullfile(qdir,fname);
        tbl = parse_tbl(dirFname); %load in summary table for each compound
        rnkPertDesc = tbl.pert_desc;
        rnkPertID = tbl.pert_id;
        tpPert = rnkPertDesc(1:1000); %top 1000 listed are positivly enriched
        tpPertID = rnkPertID(1:1000); 
        ssQ = tbl.distil_ss;

        imatch = strmatch(cmpd,tpPert); %find compound name in enrichment rank list
        n_matches(i) = length(imatch);
        
        imatchID = strmatch(cmpdID_shrt,tpPertID); %match based on pert_id
        n_matchesID(i) = length(imatchID);        
        
        sString = sprintf('%s_%sum_%s_%s',cmpd,num2str(dose1{1}),cellLine,timeP);
        sString(sString == '.') = []; %remove decimal
        sString(sString == '-') = []; %remove dash
        eString1 = sprintf('rnkStruc.%s = imatch;',sString); %save self rank events
        eString2 = sprintf('ssStruc.%s = ssQ;',sString); %save ss of queries
        eval(eString1);
        eval(eString2);
    end
    
    idiff = find(n_matches ~= n_matchesID); % where does searching by pert_desc vs pertID differ
    
%% identify which compounds had self-connection 
    iConnect = find(n_matchesID >=1);
    qConnect = pertList(iConnect);
    %qConnect = tblList(iConnect);
    qDose = cell2mat(doseList(iConnect)); 
    qnMatches = n_matchesID(iConnect)';
    
%     matchTable(:,1) = qConnect;
%     matchTable(:,2) = qDose;
%     matchTable(:,3) = qnMatches;

    %set dir for connection output
     outdir = sprintf('/xchip/cogs/projects/ASG_dose_time/cmap_queries/%s_%s_selfConnect_counts',cellLine,timeP);
     mkdir(outdir); 
     %self connection list
    [uPert,~,iupert] = unique(pertList);
    for i = 2:length(uPert)
       iCmpd = find(iupert == i);
       cmpds = pertList(iCmpd);
       doses = cell2mat(doseList(iCmpd));
       matches = n_matchesID(iCmpd)'; 
       [oDose,idose] = sort(doses); %sort doses incrementally 
       orderMatch = matches(idose); %matches in dose order
           %print match list for each compound
           outname = sprintf('%s_query_connections_isGold.txt', uPert{i}); %file name is the compound
           out = fullfile(outdir,outname); %defined above
           strings = num2cell(orderMatch);
           fid = fopen(out,'w');
           for row = 1:size(strings,1)
               fprintf(fid,'%s\n', num2str(strings{row}));
           end
           fclose(fid);
    end
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
        
        %rank plot
            outdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_24H_selfConnect_counts/viz';
            n_cmpds = length(rSigID);
            h = figure;
            hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
                         'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
                         'XLim',[1 n_cmpds],...               %#   set the x axis limit,
                         'YLim',[0 eps],...
                         'Color','none');               %#   and don't use a background color
            plot(matchList,0,'b.','MarkerSize',5);  %# Plot data set 2
            strTitle = rdb.cid{i}
            title(strrep(strTitle,'_','\_'))
            outName = sprintf('%s.jpg',lower(strTitle));
            fout = fullfile(outdir,outName);
            %control output image
                r = 150; % pixels per inch
                set(h, 'PaperUnits', 'inches', 'PaperPosition', [0 0 800 200]/r); %dim of image
                print(h,'-dpng',sprintf('-r%d',r), fout);
                %saveas(h,fout)
    end
    
    
%% get annotations for larger affogato rank
        % sort query scores
    [sMat,iMat] = sort(rdb.mat,1,'descend');
    
    n_aff = size(iMat,1); %number of perts in db
    percent = .01; %select percentage of db to examine
    n_percent = n_aff*percent; %number of compounds that represent n percent of db
    
    itmp = iMat(:,85);
    iShort = itmp(1:1000);
    sigList = rdb.rid(iShort);
    antStr = sig_info(sigList);
    
    affCmpdRnk = extractfield(antStr,'pert_desc')';
    tblCmpdRnk = tbl.pert_desc;
    affIsGld = extractfield(antStr,'is_gold')';
    
    cmpd = pertList(85);
    
    %n_is_gold = 200000; 
    n_percent = n_is_gold*percent; %number of compounds that represent n percent of is_gold
    
    %% loop through each of the compounds
    for j = 1:length(pertList)
        sScore = sMat(:,j);
        %if there was a connection, print red
        if any(j==iConnect)
            plot(sScore,'r')
            legend('connected')
            hold on
        else
            plot(sScore,'b')
            hold on
                xlabel('perturbation in affogato')
                ylabel('es score')
                titleString = sprintf('%s - %s - ASG compound connections',cellLine,timeP);
                title(titleString)
                legend('not connected')
        end
    end