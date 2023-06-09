% look at CST neural trajectories
clc
clear
close all
%% Set up
    %dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'D:\OneDrive - Northeastern University\Action Lab\01 Projects\Batista Collaboration\00 CST\01 Experiment\data\FordEarl';
        %dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
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

%%
    lambda_to_use = 3.3;
    num_dims = 8;
    start_time = -0.45;
    end_time = 0.4;

    smoothsigs = true;
    softnorm = true;
    
    % preprocess data
    td = td_preproc;
    

    if smoothsigs
        td = smoothSignals(td,struct('signals','M1_spikes','width',0.05,'calc_rate',true));
    end

    if softnorm
        td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    end

    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    % find condition independent signal for CO trials
    td_co = trimTD(td_co,{'idx_goCueTime',floor(start_time/td_co(1).bin_size)},{'idx_goCueTime',floor(end_time/td_co(1).bin_size)});
    
    
    [td_co,dpca_info] = runDPCA(td_co,'tgtDir',struct(...
        'signals','M1_spikes',...
        'num_dims',[2 6],...
        'do_plot',false,...
        'marg_names',{{'time','target'}},...
        'out_sig_prefix','M1_dpca'));
%     co_fr = cat(3,td_co.M1_spikes);
%     co_fr = mean(co_fr,3);
%     co_fr = co_fr';
%     [W, V, which_marg] = dpca( co_fr, 2, 'combinedParams', {{1}} );
%     expl_var = dpca_explainedVariance(co_fr, W, V, ...
%         'combinedParams', {{1}});
    
    % make plot of CIS activity
    timevec = start_time:td_co(1).bin_size:end_time;
    figure('defaultaxesfontsize',10)
    M1_time = cat(3,td_co.M1_dpca_time);
    indiv_dim_subplots = [1,3];
    for dim_to_plot = 1:length(indiv_dim_subplots)
        subplot(2,2,indiv_dim_subplots(dim_to_plot))
        dim_activity = squeeze(M1_time(:,dim_to_plot,:))';
        plot(timevec,dim_activity,'k','linewidth',0.5)
        hold on
        plot(timevec,mean(dim_activity),'r','linewidth',2)
        plot([0 0],get(gca,'ylim'),'g')
        set(gca,'box','off','tickdir','out','xtick',[-0.4 0 0.4])
        ylabel(sprintf('CIS dim %d',dim_to_plot))
    end
    xlabel('Time from movement onset (s)')

    subplot(2,2,[2,4])
    plot(squeeze(M1_time(:,1,:)),squeeze(M1_time(:,2,:)),'k','linewidth',0.5)
    hold on
    plot(mean(M1_time(:,1,:),3),mean(M1_time(:,2,:),3),'r','linewidth',2)
    xlabel('CIS dim 1')
    ylabel('CIS dim 2')
    set(gca,'box','off','tickdir','out','xtick',[],'ytick',[],'DataAspectRatio',[1 1 1])

    sgtitle(sprintf('%s %s CO CIS dimensions',file_query.monkey,file_query.date))

    % project CST data into CIS dimensions?