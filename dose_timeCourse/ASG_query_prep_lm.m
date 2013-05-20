%% load in infered data

    %fname = '/xchip/obelix/pod/brew/pc/ASG001_PC3_6H/by_pert_id_pert_dose/ASG001_PC3_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
    %fname = '/xchip/obelix/pod/brew/pc/ASG001_PC3_24H/by_pert_id_pert_dose/ASG001_PC3_24H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
    fname = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_6H/by_pert_id_pert_dose/ASG001_MCF7_6H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
    %fname = '/xchip/obelix/pod/brew/pc/ASG001_MCF7_24H/by_pert_id_pert_dose/ASG001_MCF7_24H_COMPZ.MODZ_SCORE_LM_n85x978.gctx';
    
    db = parse_gctx(fname);
%% set info about setup
    cellLine = 'MCF7';
    timeP = '6H';
    % define and reset gmt file names
    outdir = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/MCF7_sigs';
    outnameUP = sprintf('%s_%s_signatures_LM_up.gmt',cellLine,timeP);
    outnameDN = sprintf('%s_%s_signatures_LM_dn.gmt',cellLine,timeP);
    
    %delete contents of existing UP file
    out = fullfile(outdir,outnameUP);
    fid = fopen(out,'w');
    fclose(fid);
    
    %delete contents of existing DN file
    out = fullfile(outdir,outnameDN);
    fid = fopen(out,'w');
    fclose(fid);
    
%% loop though all signatures in plate - save signature
pertList = db.cdesc(:,db.cdict('pert_desc')); %identify DMSO
for i = 1:length(pertList);
% for i = 1:4;

    iCmpd = i; %alvespimycinv at 10um
    cmpd = db.cdesc(iCmpd,db.cdict('pert_desc'));
    dose1 = db.cdesc(iCmpd,db.cdict('pert_dose'));
    profile = db.mat(:,iCmpd);
    [zSort, ia] = sort(profile,'descend');
    srtProbe = db.rid(ia);
    
    %isolate edges
    n_edge = 50;
    n_probes = length(srtProbe);
    probe_up = srtProbe(1:n_edge);
    probe_dn = srtProbe((n_probes-n_edge+1):n_probes);
    
    % create gmt files - NOTE: THIS CODE WRITES AN EXTRA SPACE WHICH MUST
    % BE RMEOVED FROM TOP OF THE GMT FILE
        %append signature to UP gmt file
        geneSetName = sprintf('%s_%sum_%s_%s',cmpd{1},num2str(dose1{1}),cellLine,timeP);
        geneSetDesc = geneSetName;
        out = fullfile(outdir,outnameUP); %defined above
        strings = probe_up;
        fid = fopen(out,'a');
        fprintf(fid, '\n%s\t%s\t', geneSetName,geneSetDesc); %start the gmt file by printing set name/ desc
        for row = 1:size(strings,1)
            fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
            fprintf(fid, '%s\t', strings{row,end});
        end
        fclose(fid);    
    
        %append signature to DN gmt file
        geneSetName = sprintf('%s_%sum_%s_%s',cmpd{1},num2str(dose1{1}),cellLine,timeP);
        geneSetDesc = geneSetName;
        out = fullfile(outdir,outnameDN); %defined above
        strings = probe_dn;
        fid = fopen(out,'a');
        fprintf(fid, '\n%s\t%s\t', geneSetName,geneSetDesc);
        for row = 1:size(strings,1)
            fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
            fprintf(fid, '%s\t', strings{row,end});
        end
        fclose(fid);
    
    
%         %write individual signature files
%         outnameUP = sprintf('%s_%sum_%s_%s_up.txt',cmpd{1},num2str(dose1{1}),cellLine,timeP);
%         %print UP list to a file
%         out = fullfile(outdir,outnameUP); %defined above
%         strings = probe_up;
%         fid = fopen(out,'w');
%         for row = 1:size(strings,1)
%             fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
%             fprintf(fid, '%s\n', strings{row,end});
%         end
%         fclose(fid);
%     
%         %print DN list to a file
%         outnameDN = sprintf('%s_%sum_%s_%s_dn.txt',cmpd{1},num2str(dose1{1}),cellLine,timeP);
%         out = fullfile(outdir,outnameDN); %defined above
%         strings = probe_dn;
%         fid = fopen(out,'w');
%         for row = 1:size(strings,1)
%             fprintf(fid, repmat('%s\t',1,size(strings,2)-1), strings{row,1:end-1});
%             fprintf(fid, '%s\n', strings{row,end});
%         end
%         fclose(fid);
end
