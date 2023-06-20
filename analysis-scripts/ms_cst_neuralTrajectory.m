    % look at CST neural trajectories
    clc
    clear
    close all
    set(groot,'defaultAxesFontSize',13)

    %% Set up

    dataroot = '/Users/Shared/Previously Relocated Items/Security/00 Postdoc at Northeastern/01 Projects/Batista Collaboration/00 CST /cst-gainlag_Raeed/All data';
    if ispc
        dataroot = 'D:\OneDrive - Northeastern University\Action Lab\01 Projects\Batista Collaboration\00 CST\01 Experiment\data\FordEarl';
    end

    file_info = dir(fullfile(dataroot,'*COCST*'));
    filenames = horzcat({file_info.name})';
    fileDates = {'2019/07/02' , '2019/07/16' , '2018/06/18' , '2018/06/26' , '2018/06/27'};
    monkeynames = {'Earl','Earl','Ford','Ford','Ford'};

    %% Select a file
    
    ID = 1; % monkey ID: Earl:1,2 Ford:3,4,5
    
    file_query = struct(...
        'monkey','Ford',...
        'date','20180627');
    %td_preproc = load_clean_cst_data(fullfile(dataroot,sprintf('%s_%s_COCST_TD.mat',file_query.monkey,file_query.date)));
    td_preproc = load_clean_cst_data(fullfile(dataroot,sprintf('%s',filenames{ID})));

    % Make sure we have CST trials
    assert(~isempty(td_preproc),sprintf('Incomplete dataset for file %s %s\n', file_query.monkey,file_query.date))
    assert(~isempty(td_preproc(1).M1_unit_guide),sprintf('Skipping file %s %s because no spike data...\n',file_query.monkey,file_query.date))

    %% Smoothing the data
    num_dims = 8;
    softnorm = true;
    smoothsigs = true;
    td = td_preproc;

    if smoothsigs
        td = smoothSignals(td,struct('signals','M1_spikes','width',0.050,'calc_rate',true));
    end

    % trim TD to only center hold portion
    startSample = 0;
    endSample   = 400;
    td = trimTD(td,{'idx_goCueTime',startSample},{'idx_goCueTime',endSample});

    % smooth data
    if softnorm
        td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    end


    %% Find trial index for control strategy
    Labels = load('Labels'); Labels = Labels.Labels;
    
    iip = contains(Labels.monkey,monkeynames{ID}) & ...
         contains(Labels.date,fileDates{ID}) & ...
         contains(Labels.ControlPolicy,'Position');
    PosTrlId =  Labels.trial_id(iip);
    
    
    iiv = contains(Labels.monkey,monkeynames{ID}) & ...
         contains(Labels.date,fileDates{ID}) & ...
         contains(Labels.ControlPolicy,'Velocity');
    VelTrlId =  Labels.trial_id(iiv);
    
        
    iih = contains(Labels.monkey,monkeynames{ID}) & ...
         contains(Labels.date,fileDates{ID}) & ...
         contains(Labels.ControlPolicy,'Hybrid');
    HybTrlId =  Labels.trial_id(iih);
    

    
    %%  Dimensionality reduction
    % td = dimReduce(td,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));

    % get average for each bin
    td_binned = binTD(td,'average');
    %td_binned = binTD(td,10);

    td_binned = dimReduce(td_binned,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));

    % split data
    [~,td_co]  = getTDidx(td_binned,'task','CO');
    [~,td_cst] = getTDidx(td_binned,'task','CST');
    %[~,td_co] = getTDidx(td_binned,'task','CST','lambda',4.1);

    
    CstTrlId = [td_cst.trial_id]';
    
    
    %% CMD Scaling / PCA
    
    M1 = cat(1,td_cst.M1_spikes);
    Y = pdist(M1);
    Z = squareform(Y);
    [Yc , ic] = cmdscale(Z,5);

    
    
    figure(1)
    clf
    sb1=2;
    sb2=4;
    GroupCl = [.8,.5,.0 ; .0,.7,.7 ; .7,.7,.7];
    subplot(sb1,sb2,1)
    ip = ismember(CstTrlId , PosTrlId);
    iv = ismember(CstTrlId , VelTrlId);
    ih = ismember(CstTrlId , HybTrlId);
    plot3(Yc(ip,3),Yc(ip,2),Yc(ip,1),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot3(Yc(iv,3),Yc(iv,2),Yc(iv,1),'.','color',GroupCl(2,:),'markersize',16);
    plot3(Yc(ih,3),Yc(ih,2),Yc(ih,1),'.','color',GroupCl(3,:),'markersize',16);
    grid
    xlabel('Config 3')
    ylabel('Config 2')
    zlabel('Config 1')
    
    subplot(sb1,sb2,2)
    plot(Yc(ip,2),Yc(ip,1),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot(Yc(iv,2),Yc(iv,1),'.','color',GroupCl(2,:),'markersize',16);
    plot(Yc(ih,2),Yc(ih,1),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('Config 2')
    ylabel('Config 1')
    legend('Pos','Vel','None')
    box off

    subplot(sb1,sb2,3)
    plot(Yc(ip,3),Yc(ip,1),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot(Yc(iv,3),Yc(iv,1),'.','color',GroupCl(2,:),'markersize',16);
    plot(Yc(ih,3),Yc(ih,1),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('Config 3')
    ylabel('Config 1')
    legend('Pos','Vel','None')
    box off
    
    subplot(sb1,sb2,4)
    plot(Yc(ip,3),Yc(ip,2),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot(Yc(iv,3),Yc(iv,2),'.','color',GroupCl(2,:),'markersize',16);
    plot(Yc(ih,3),Yc(ih,2),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('Config 3')
    ylabel('Config 2')
    legend('Pos','Vel','None')
    box off
    
    % PCA analysis
    
    M1_pca = cat(1,td_cst.M1_pca);
    
    subplot(sb1,sb2,sb2+1)
    plot3(M1_pca(ip,3),M1_pca(ip,2),M1_pca(ip,1),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot3(M1_pca(iv,3),M1_pca(iv,2),M1_pca(iv,1),'.','color',GroupCl(2,:),'markersize',16);
    plot3(M1_pca(ih,3),M1_pca(ih,2),M1_pca(ih,1),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('PC3')
    ylabel('PC2')
    zlabel('PC1')
    grid

    subplot(sb1,sb2,sb2+2)
    plot(M1_pca(ip,2),M1_pca(ip,1),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot(M1_pca(iv,2),M1_pca(iv,1),'.','color',GroupCl(2,:),'markersize',16);
    plot(M1_pca(ih,2),M1_pca(ih,1),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('PC2')
    ylabel('PC1')
    legend('Pos','Vel','None')
    box off
    
    subplot(sb1,sb2,sb2+3)
    plot(M1_pca(ip,3),M1_pca(ip,1),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot(M1_pca(iv,3),M1_pca(iv,1),'.','color',GroupCl(2,:),'markersize',16);
    plot(M1_pca(ih,3),M1_pca(ih,1),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('PC3')
    ylabel('PC1')
    legend('Pos','Vel','None')
    box off
    
    subplot(sb1,sb2,sb2+4)
    plot(M1_pca(ip,3),M1_pca(ip,2),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot(M1_pca(iv,3),M1_pca(iv,2),'.','color',GroupCl(2,:),'markersize',16);
    plot(M1_pca(ih,3),M1_pca(ih,2),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('PC3')
    ylabel('PC2')
    legend('Pos','Vel','None')
    box off
    
    return
    
    %%

    % make plot of hold time activity
    figure(2)
    clf
    pax1 = 1; % which PCA on axis 1?
    pax2 = 2;
    pax3 = 3;

    subplot(1,2,1)
    scatter3(...
        get_vars(td_co,{'M1_pca',pax1}),...
        get_vars(td_co,{'M1_pca',pax2}),...
        get_vars(td_co,{'M1_pca',pax3}),...
        [],'r','filled');
    hold on
    scatter3(...
        get_vars(td_cst,{'M1_pca',pax1}),...
        get_vars(td_cst,{'M1_pca',pax2}),...
        get_vars(td_cst,{'M1_pca',pax3}),...
        [],cat(1,td_cst.lambda),'filled');
    colormap(viridis);
    set(gca,'box','off','tickdir','out')
    axis equal
    xlabel('PC1')
    ylabel('PC2')
    zlabel('PC3')
    title('M1 PCA - CO in red, CST in viridis')

    subplot(1,2,2)
    scatter3(...
        get_vars(td_co,{'hand_pos',1}),...
        get_vars(td_co,{'hand_pos',2}),...
        get_vars(td_co,{'hand_pos',3}),...
        [],'r','filled');
    %     scatter3(...
    %         get_vars(td_co,{'rel_hand_pos',1}),...
    %         get_vars(td_co,{'rel_hand_pos',2}),...
    %         get_vars(td_co,{'rel_hand_pos',3}),...
    %         [],'r','filled');
    hold on
    scatter3(...
        get_vars(td_cst,{'hand_pos',1}),...
        get_vars(td_cst,{'hand_pos',2}),...
        get_vars(td_cst,{'hand_pos',3}),...
        [],cat(1,td_cst.lambda),'filled');
    %     scatter3(...
    %         get_vars(td_cst,{'rel_hand_pos',1}),...
    %         get_vars(td_cst,{'rel_hand_pos',2}),...
    %         get_vars(td_cst,{'rel_hand_pos',3}),...
    %         [],cat(1,td_cst.lambda),'filled');
    set(gca,'box','off','tickdir','out')
    axis equal
    xlabel('Hand X position')
    ylabel('Hand Y position')
    zlabel('Hand Z position')
    title('Hand position - CO in red, CST in viridis')

    if smoothsigs
        smooth_str = '50ms';
    else
        smooth_str = 'no';
    end
    sgtitle(sprintf('%s %s hold period activity (%s smooth)',file_query.monkey,file_query.date,smooth_str))

    %% plot out hold time behavior
    co_idx = getTDidx(td,'task','CO');
    cst_idx = getTDidx(td,'task','CST');
    co_hand_pos = cat(1,td(co_idx).hand_pos);
    cst_hand_pos = cat(1,td(cst_idx).hand_pos);
    %     co_hand_pos = cat(1,td(co_idx).rel_hand_pos);-+
    %     cst_hand_pos = cat(1,td(cst_idx).rel_hand_pos);
    figure('defaultaxesfontsize',18)
    plot(cst_hand_pos(:,1),cst_hand_pos(:,2))
    hold on
    plot(co_hand_pos(:,1),co_hand_pos(:,2))
    axis equal
    set(gca,'box','off','tickdir','out')
    legend('CST hold','CO hold')
    xlabel('Hand x-position')
    ylabel('Hand y-position')
