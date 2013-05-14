%% load sig ID file
%     fname = '/xchip/cogs/projects/dos/mongo_queries/DOS_sig_ids.txt';
%     sigList = importdata(fname);

%     fname = '/xchip/cogs/projects/dos/mongo_queries/DOS_sig_ids+isGld.txt';
%     sigList = importdata(fname,' ');
    
        %% loop through each sig ID
%         for i = 1:length(sigList);
%             annt = sig_info(sigList{1});
%             goldList(i) = annt.is_gold;
%             cellList{i} = annt.cell_id;
%         end
%     

%% load data and define columns
    fname = '/xchip/cogs/projects/dos/mongo_queries/DOS_sig_info.txt';
    sigList = importdata(fname,' ');
    
    sig_id = sigList.textdata(:,1);
    isGld = sigList.textdata(:,2);
    isGld = cell2mat(isGld); %convert to mat
    isGld = str2num(isGld); %convert to num
    cellID = sigList.textdata(:,3);
    pertID = sigList.textdata(:,4);
    ss = sigList.data;

%% evaluate which experiments are is gold    
    ig = find(isGld == 1); %get index of isGld
    n_isGld = sum(isGld(ig));
    gldCell = cellID(ig); %what are the cell lines for each gold instance
    gldSS = ss(ig);
    gldPert = pertID(ig);
    gldSigId = sig_id(ig);
%% find perturbation that were gold in more than once in the db
    [unGldPert,~,igldPert] = unique(gldPert); % find the unique DOS pertIDs that ar is_gold
    [gldPerCmpd,~]=hist(igldPert,unique(igldPert)); %count how many times a pertID occurs as is_gold
    
        [n1, xout1] = hist(gldPerCmpd,1:1:6,100);
        bar(xout1,n1,'b'); grid;
            xlabel('number of occurences')
            ylabel('freq')
            title('unique DOS compounds as isGold')      
    
    %does average ss increase with with multiple occurences
    avSSlist = nan(1,6);
    for j = 3:length(xout1)
        j =6
        n_occ = j; %occurrence cuttoff        
        iOccN = find(gldPerCmpd == n_occ); %find occurences in unique list
        OccName = unGldPert(iOccN); %Pert Names of multiple occurences
                
        % loop through each set of compounds
        clear OssList
        for i = 1:length(OccName)
            cmpd = OccName(i); %which item in the list
            iOall = strmatch(cmpd,pertID);
            iOgld = strmatch(cmpd,gldPert);

            OssList(i) = mean(gldSS(iOgld)); %averae ss
            Ocell = gldCell(iOgld); %cell lines of occurences
        end

        avSSlist(j) = mean(OssList);
    end
%% get sig Ids of specific compounds
    one_cmpd = OccName;
    imtch = strmatch(one_cmpd, pertID);
    cmpd_sigID = sig_id(imtch);
    
%% what are the perturbation parameters: cell lines, dose, time point
    [unPert,~,iPert] = unique(pertID); % find the unique perts
    [cntPerCmpd,~]=hist(iPert,unique(iPert)); %count how many times a pertID occurs in ds
    
        [n1, xout1] = hist(cntPerCmpd,1:1:7);
        bar(xout1,n1,'b'); grid;
            xlabel('number of occurences')
            ylabel('freq')
            title('unique DOS compounds in affogato')     
%% cell lines
    [unDosCell,~,iUnDCell] = unique(cellID); %unique DOS cell lines
    [unGldCell,~,iUnCell] = unique(gldCell); %unique cell lines with isGld
   
%% plot sig strength 
        [n1, xout1] = hist(ss,0:.4:20,100);
        bar(xout1,n1,'r'); grid; hold on
        [n2, xout2] = hist(gldSS,0:.4:20,100);
        bar(xout2,n2,'g');
        hold off
            xlabel('sig strength')
            ylabel('freq')
            title('Signature strength for DOS compounds')
            legend('all DOS','only isGold')

%% plot cell line distribution
    hist(iUnDCell,6)
        [n1, xout1] = hist(iUnDCell,1:1:6,100);
        bar(xout1,n1,'r'); grid; hold on
        [n2, xout2] = hist(iUnCell,1:1:6,100);
        bar(xout2,n2,'g');
        hold off
            xlabel('cell line')
            ylabel('freq')
            title('DOS cell line distribution')
            legend('all DOS','isGold')
            set(gca,'xticklabel',unGldCell)
    
    