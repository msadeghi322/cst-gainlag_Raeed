% look at CST neural trajectories
clc
clear
close all

%% Set up
    if ispc
        dataroot = 'D:\OneDrive - Northeastern University\Action Lab\01 Projects\Batista Collaboration\00 CST\01 Experiment\data\FordEarl';
    end

    file_info = dir(fullfile(dataroot,'*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% Select a file
    file_query = struct(...
        'monkey','Ford',...
        'date','20180627');
    td_preproc = load_clean_cst_data(fullfile(dataroot,sprintf('%s_%s_COCST_TD.mat',file_query.monkey,file_query.date)));
    
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
    startSample = -400;
    endSample   = 0;
    td = trimTD(td,{'idx_goCueTime',startSample},{'idx_goCueTime',endSample});

    % smooth data
    if softnorm
        td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    end

    
    [~,td_cst] = getTDidx(td,'task','CST');
    M1 = cat(3,td_cst.M1_spikes);
    figure(1)
    clf
    sb1=2;
    sb2=2;
    subplot(sb1,sb2,1)
    y = squeeze(M1(:,20,:));
    t = (startSample:endSample)*td_cst(1).bin_size;
    plot(t,y,'color',.7*ones(1,3));
    
    
%%  Dimensionality reduction
    % td = dimReduce(td,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));

    % get average for each bin
    td_binned = binTD(td,'average');
    %td_binned = binTD(td,10);

    td_binned = dimReduce(td_binned,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
    
    % split data
    [~,td_co] = getTDidx(td_binned,'task','CO');
    [~,td_cst] = getTDidx(td_binned,'task','CST');
    %[~,td_co] = getTDidx(td_binned,'task','CST','lambda',4.1);
    
    
    
    
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
