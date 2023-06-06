% look at CST neural trajectories

%% Set up
    %dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'D:\OneDrive - Northeastern University\Action Lab\01 Projects\Batista Collaboration\00 CST\01 Experiment\data\FordEarl';
    end

    file_info = dir(fullfile(dataroot,'*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% Loop through files
    file_query = struct(...
        'monkey','Earl',...
        'date','20190716');
    td_preproc = load_clean_cst_data(fullfile(dataroot,sprintf('%s_%s_COCST_TD.mat',file_query.monkey,file_query.date)));
    
    % Make sure we have CST trials
    assert(~isempty(td_preproc),sprintf('Incomplete dataset for file %s %s\n', file_query.monkey,file_query.date))
    assert(~isempty(td_preproc(1).M1_unit_guide),sprintf('Skipping file %s %s because no spike data...\n',file_query.monkey,file_query.date))
    
%% 
    % options
    softnorm = true;
    analysis_epoch = 'move';
    smoothsigs = false;

    td = td_preproc;
    
    % smooth data
    if smoothsigs
        td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
    end
    
    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    lambdas = sort(unique(cat(1,td_cst.lambda)));
    
    % trim TD to some arbitrary window after go cue
    if strcmpi(analysis_epoch,'move')
        td_co = trimTD(td_co,{'idx_goCueTime',150},{'idx_goCueTime',500});
        td_cst = trimTD(td_cst,{'idx_cstStartTime',150},{'idx_cstStartTime',5000});
    elseif strcmpi(analysis_epoch,'hold')
        td_co = trimTD(td_co,{'idx_goCueTime',-450},{'idx_goCueTime',0});
        td_cst = trimTD(td_cst,{'idx_goCueTime',-450},{'idx_goCueTime',0});
    elseif strcmpi(analysis_epoch,'full')
        td_co = trimTD(td_co,{'idx_goCueTime',-400},{'idx_goCueTime',500});
        td_cst = trimTD(td_cst,{'idx_goCueTime',-400},{'idx_goCueTime',5000});
    end

    td_co = binTD(td_co,50);
    td_cst = binTD(td_cst,50);

    % soft normalize
    if softnorm
        td_co = softNormalize(td_co,struct('signals','M1_spikes'));
        td_cst = softNormalize(td_cst,struct('signals','M1_spikes'));
    end
    
    % Get variances
    total_co_var = sum(var(getSig(td_co,'M1_spikes')));
    total_co_pr = calc_participation_ratio(getSig(td_co,'M1_spikes'));
    total_co_pa_dim = parallel_analysis_dimensionality(getSig(td_co,'M1_spikes'),struct('do_plot','true'));
    total_cst_var = sum(var(getSig(td_cst,'M1_spikes')));
    total_cst_pr = calc_participation_ratio(getSig(td_cst,'M1_spikes'));
    total_cst_pa_dim = parallel_analysis_dimensionality(getSig(td_cst,'M1_spikes'),struct('do_plot','true'));

    % for each lambda...
    lambda_var = zeros(length(lambdas),1);
    lambda_pr = zeros(length(lambdas),1);
    lambda_pa_dim = zeros(length(lambdas),1);
    for lambdanum = 1:length(lambdas)
        [~,td_lambda] = getTDidx(td_cst,'lambda',lambdas(lambdanum));
        lambda_var(lambdanum) = sum(var(getSig(td_lambda,'M1_spikes')));
        lambda_pr(lambdanum) = calc_participation_ratio(getSig(td_lambda,'M1_spikes'));
        lambda_pa_dim(lambdanum) = parallel_analysis_dimensionality(getSig(td_lambda,'M1_spikes'));
    end

    % for random sets of CST trials with the same total amount of data as CO
    num_shuffles = 45;
    num_cst_trials = floor(length(td_co)/10);
    shuffle_var = zeros(num_shuffles,1);
    shuffle_pr = zeros(num_shuffles,1);
    shuffle_pa_dim = zeros(num_shuffles,1);
    for shufflenum = 1:num_shuffles
        [~,td_shuffle] = getTDidx(td_cst,'rand',num_cst_trials);
        shuffle_var(shufflenum) = sum(var(getSig(td_shuffle,'M1_spikes')));
        shuffle_pr(shufflenum) = calc_participation_ratio(getSig(td_shuffle,'M1_spikes'));
        shuffle_pa_dim(shufflenum) = parallel_analysis_dimensionality(getSig(td_shuffle,'M1_spikes'));
    end
    
    %% plot out summaries
    figure('defaultaxesfontsize',18)
    colors = linspecer(2);
    cst_color = colors(1,:);
    co_color = colors(2,:);

    % plot variancess
    subplot(3,1,1)
    plot(lambdas,lambda_var,'-','linewidth',2,'color',cst_color)
    hold on
    plot(lambdas,repmat(total_co_var,size(lambdas)),'--r','linewidth',2,'color',co_color)
    plot(lambdas,repmat(total_cst_var,size(lambdas)),'--b','linewidth',2,'color',cst_color)
    plot(lambdas,repmat(mean(shuffle_var),size(lambdas)),'--','linewidth',2,'color',cst_color*0.5)
    title('Neural variance')
    set(gca,'box','off','tickdir','out','xlim',[min(lambdas) max(lambdas)],'xtick',[])
    
    % plot PRs
    subplot(3,1,2)
    plot(lambdas,lambda_pr,'-','linewidth',2,'color',cst_color)
    hold on
    plot(lambdas,repmat(total_co_pr,size(lambdas)),'--','linewidth',2,'color',co_color)
    plot(lambdas,repmat(total_cst_pr,size(lambdas)),'--','linewidth',2,'color',cst_color)
    plot(lambdas,repmat(mean(shuffle_pr),size(lambdas)),'--','linewidth',2,'color',cst_color*0.5)
    title('Dimensionality (via PR)')
    set(gca,'box','off','tickdir','out','xlim',[min(lambdas) max(lambdas)],'xtick',[])
    
    % plot PA dims
    subplot(3,1,3)
    plot(lambdas,lambda_pa_dim,'-','linewidth',2,'color',cst_color)
    hold on
    plot(lambdas,repmat(total_co_pa_dim,size(lambdas)),'--','linewidth',2,'color',co_color)
    plot(lambdas,repmat(total_cst_pa_dim,size(lambdas)),'--','linewidth',2,'color',cst_color)
    plot(lambdas,repmat(mean(shuffle_pa_dim),size(lambdas)),'--','linewidth',2,'color',cst_color*0.5)
    xlabel('\lambda')
    title('Dimensionality (via parallel ananlysis)')
    set(gca,'box','off','tickdir','out','xlim',[min(lambdas) max(lambdas)])

    sgtitle(sprintf('%s %s %s time',file_query.monkey,file_query.date,analysis_epoch))
