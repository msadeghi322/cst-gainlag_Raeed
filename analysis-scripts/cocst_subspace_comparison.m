clc 
clear
close all

%% setup
    %dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'D:\OneDrive - Northeastern University\Action Lab\01 Projects\Batista Collaboration\00 CST\01 Experiment\data\FordEarl';
    end
    
%% Load a file
    file_query = struct(...
        'monkey','Earl',...
        'date','20190716');
    td_preproc = load_clean_cst_data(fullfile(dataroot,sprintf('%s_%s_COCST_TD.mat',file_query.monkey,file_query.date)));
    
    % Make sure we have CST trials
    assert(~isempty(td_preproc),sprintf('Incomplete dataset for file %s %s\n', file_query.monkey,file_query.date))
    assert(~isempty(td_preproc(1).M1_unit_guide),sprintf('Skipping file %s %s because no spike data...\n',file_query.monkey,file_query.date))

%% Preprocess td some more...
    % options
    analysis_epoch = 'move';
    smoothsigs = true;
    num_dims = 20;
    equalize_samples = false;
    num_samples = 1500;
    avg_co_data = false;
%     assignParams(who,params);

    td = td_preproc;
    
    % smooth data
    if smoothsigs
        td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true,'field_extra','_smooth'));
    end
    
    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    % trim TD to some arbitrary window after go cue
    if strcmpi(analysis_epoch,'move')
        td_co = trimTD(td_co,struct('idx_start',{{'idx_goCueTime',150}},'idx_end',{{'idx_goCueTime',500}},'remove_short',true));
        td_cst = trimTD(td_cst,struct('idx_start',{{'idx_cstStartTime',150}},'idx_end',{{'idx_cstStartTime',5000}},'remove_short',true));
    elseif strcmpi(analysis_epoch,'hold')
        td_co = trimTD(td_co,struct('idx_start',{{'idx_goCueTime',-450}},'idx_end',{{'idx_goCueTime',0}},'remove_short',true));
        td_cst = trimTD(td_cst,struct('idx_start',{{'idx_goCueTime',-450}},'idx_end',{{'idx_goCueTime',0}},'remove_short',true));
    end

    td_co = binTD(td_co,50);
    td_cst = binTD(td_cst,50);

    [td_co,pca_co] = dimReduce(td_co,struct(...
        'signals','M1_spikes_smooth',...
        'algorithm','pca',...
        'num_dims',num_dims));
    [td_cst,pca_cst] = dimReduce(td_cst,struct(...
        'signals','M1_spikes_smooth',...
        'algorithm','pca',...
        'num_dims',num_dims));

    [td_co.M1_pca_co] = deal(td_co.M1_smooth_pca);
    [td_cst.M1_pca_cst] = deal(td_cst.M1_smooth_pca);

    td_co = dimReduce(td_co,pca_cst);
    [td_co.M1_pca_cst] = deal(td_co.M1_smooth_pca);
    td_cst = dimReduce(td_cst,pca_co);
    [td_cst.M1_pca_co] = deal(td_cst.M1_smooth_pca);

    td_co = rmfield(td_co,'M1_smooth_pca');
    td_cst = rmfield(td_cst,'M1_smooth_pca');

    total_co_var = sum(var(cat(1,td_co.M1_pca_co)));
    total_cst_var = sum(var(cat(1,td_cst.M1_pca_cst)));

    percent_co_var_in_co_space = var(cat(1,td_co.M1_pca_co))/total_co_var;
    percent_co_var_in_cst_space = var(cat(1,td_co.M1_pca_cst))/total_co_var;
    percent_cst_var_in_cst_space = var(cat(1,td_cst.M1_pca_cst))/total_cst_var;
    percent_cst_var_in_co_space = var(cat(1,td_cst.M1_pca_co))/total_cst_var;

    figure
    ax(1) = subplot(1,2,1);
    bar([percent_cst_var_in_co_space;percent_co_var_in_co_space]')
    title('CO subspace')
    legend('CST data','CO data')
    set(gca,'box','off','tickdir','out')
    ax(2) = subplot(1,2,2);
    bar([percent_cst_var_in_cst_space;percent_co_var_in_cst_space]')
    title('CST subspace')
    set(gca,'box','off','tickdir','out')
    linkaxes(ax,'y')

    