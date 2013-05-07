addpath('/xchip/cogs/hogstrom/scripts/FDR') %add scripts for FDR correction
addpath('/xchip/cogs/hogstrom/scripts/outside_code/resampling_statistical_toolkit/statistics') %second script for fdr correction
%% load in query result
%--> GBP = brewed by pert_id
% dmsoESfile = '/xchip/cogs/hogstrom/analysis/DMSO_nullDist/apr02/my_analysis.query_tool.2013040214375649/result_WTES.COMBINED_n558x398050.gctx';
% --> GBW = brewed by RNA well
dmsoESfile = '/xchip/cogs/hogstrom/analysis/DMSO_nullDist/apr04/my_analysis.query_tool.2013040415074749/result_WTES.COMBINED_n1000x398050.gctx';
dmsoRslt = parse_gctx(dmsoESfile);


% query results from 20 poscons
% posconQueryFile = '/xchip/cogs/hogstrom/analysis/DMSO_nullDist/apr03/my_analysis.query_tool.2013040310092991/result_WTES.COMBINED_n20x398050.gctx';
% poscon = parse_gctx(posconQueryFile);

% query results from complete set of poscons
posconQueryFile = '/xchip/cogs/hogstrom/analysis/DMSO_nullDist/apr05/my_analysis.query_tool.2013040516234391/result_WTES.COMBINED_n372x398050.gctx';
poscon = parse_gctx(posconQueryFile);

%% remove dmso self-connections
nDmsos = length(dmsoRslt.mat(1,:));

% ind1s = (dmsoRslt.mat == 1);
% ind1s = find(dmsoRslt.mat == 1);
for i = 1:nDmsos
    ind1s = find(dmsoRslt.mat(:,i) >= .99);
    if length(ind1s) >= 1
        annt1 = dmsoRslt.rid(ind1s);
        annt2 = dmsoRslt.cid(i);
        equalList(i) = isequal(annt1,annt2);
    end
    %set self connections to zero
    dmsoRslt.mat(ind1s,i) = 0;
    ind1s = [];
end
%% one tailed ES test
% nDmsos = length(dmsoRslt.mat(1,:));
% esCol = poscon.mat(:,1);
% % for i in 
% for i = 1:nDmsos
%     dmsoCol = dmsoRslt.mat(:,i);
%     %grtr = esCol > dmsoCol;
%     grtr = dmsoCol > esCol;
%     grtrMtrx(:,i) = grtr;
% end
% nGrtr = sum(grtrMtrx,2);
% pval = nGrtr/nDmsos;
% [sortP, iSortP] = sort(pval);

%% two tailed ES test
% nPoscon = length(poscon.mat(1,:));
% nDmsos = length(dmsoRslt.mat(1,:));
% for j = 1:nPoscon
%     esCol = poscon.mat(:,j); 
%     % for i in 
%     for i = 1:nDmsos
%         dmsoCol = dmsoRslt.mat(:,i);
%         %grtr = esCol > dmsoCol;
%         grtr = dmsoCol > esCol;
%         grtrMtrx(:,i) = grtr;
%         lss = dmsoCol < esCol;
%         lssMtrx(:,i) = lss;
%     end
%     nGrtr = sum(grtrMtrx,2);
%     nLss = sum(lssMtrx,2);
%     twoTailcount = min(nGrtr,nLss);
%     pval = twoTailcount/(nDmsos*2);
%     pvalMtrx(:,j) = pval;
%     [sortP, iSortP] = sort(pval);
%     esBySig = esCol(iSortP);
% end

%% two tailed ES test - assume uniform distribution around zero
nPoscon = length(poscon.mat(1,:));
nDmsos = length(dmsoRslt.mat(1,:));
for j = 1:nPoscon
    esCol = poscon.mat(:,j); 
    % for i in 
    for i = 1:nDmsos
        dmsoCol = dmsoRslt.mat(:,i);
        grtr = abs(dmsoCol) > abs(esCol);
        grtrMtrx(:,i) = grtr;
    end
    pval = (1+sum(grtrMtrx,2))/(nDmsos);
    pvalMtrx(:,j) = pval;
    %[sortP, iSortP] = sort(pval);
    %esBySig = esCol(iSortP);
end

%% run FDR correction
q = .1;
% [pID,pN] = FDR(pval,q);
% [p_fdr, p_masked] = fdr( pval, q);

pvalThreshMtrx = zeros(size(pvalMtrx,1),size(pvalMtrx,2));
for iposSig = 1:size(pvalMtrx,2)
    [pID,pN] = FDR(pvalMtrx(:,iposSig),q);
    if isempty(pID) == 1
        pID = 0;
    end
    FDRlist(iposSig) = pID;
    passCount = sum(pvalMtrx(:,iposSig) < pID);
    passCountList(iposSig) = passCount;
%     pvalThreshMtrx(:,pvalMtrx(:,iposSig) < pID) = pvalMtrx(pvalMtrx(:,iposSig) < pID,iposSig);
end

iHits = find(passCountList >1);
poscon.cid(iHits);

%% 
hist(passCountList,20)
xlabel('number of significant connections')
ylabel('freq')
title('372 poscon cmap queries')

%% what is the relationship between ss and average connection
%what are the es scores that pass the fdr threshold


%% find NaNs and zeros in origonal matrix
% find(isnan(esCol))
