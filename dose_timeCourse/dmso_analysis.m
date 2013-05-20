%% load in DMSO instances from the db
%     fname = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/cmpd_sigIDs/DMSO_byID.grp';
    fname = '/xchip/cogs/projects/ASG_dose_time/cmap_queries/cmpd_sigIDs/DMSO_byID_short.grp';
    gctname = '/xchip/cogs/data/build/affogato/affogato_r1_score_n398050x22268.gctx';
    
%% load z-score data into memory 
    db = parse_gctx(gctname,'cid',fname);

    %get annotations
    annot = sig_info(fname);
    
    affCmpd = extractfield(annot,'pert_desc')';
    affIsGld = extractfield(annot,'is_gold')';
    affID = extractfield(annot,'pert_id')';
    affSS = extractfield(annot,'distil_ss')';
    
%% qq plots
figure(1)
subplot(2,1,1)
qqplot(db.mat(:,1))
subplot(2,1,2)
hist(db.mat(:,1),50)

figure(2)
subplot(2,1,1)
qqplot(db.mat(:,7))
subplot(2,1,2)
hist(db.mat(:,7),50)

figure(3)
subplot(2,1,1)
qqplot(db.mat(:,100))
subplot(2,1,2)
hist(db.mat(:,100),50)

figure(4)
subplot(2,1,1)
qqplot(db.mat(:,25))
subplot(2,1,2)
hist(db.mat(:,25),50)

figure(5)
subplot(2,1,1)
qqplot(db.mat(:,55))
subplot(2,1,2)
hist(db.mat(:,55),50)
