clear
clc
close all
set(0,'DefaultFigureWindowStyle','docked')
%% Set up
    if ispc
        dataroot = 'D:\OneDrive - Northeastern University\Action Lab\01 Projects\Batista Collaboration\00 CST\01 Experiment\data\FordEarl';
       % dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
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
    
%% Smooth and extract td_cst
    td = td_preproc;
    softnorm = true;
    % smooth data
%     td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
    % softnorm
    if softnorm
        td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    end
    
%% Look for neural engagement signal in center hold period
    [~,td_cst] = getTDidx(td,'task','CST');
    td_hold = trimTD(td_cst,struct(...
        'idx_start',{{'idx_goCueTime',-400}},...
        'idx_end',{{'idx_goCueTime',0}},...
        'remove_short',true));
    td_hold = binTD(td_hold,'average');
    
    [td_hold,hold_pca_info] = dimReduce(td_hold,struct('algorithm','pca','signals','M1_spikes'));
    
    % Hold activity PC 1 plotted with indications of lambda changes and failures
    lambda_changes = find(diff(cat(1,td_hold.lambda))~=0)';
    failure_idx = getTDidx(td_hold,'result','F');
    line_colors = linspecer(2);
    figure('defaultaxesfontsize',10)
    ax(1) = subplot(4,1,1);
    imagesc(cat(1,td_hold.lambda)')
    hold on
    plot(repmat(failure_idx,2,1)+0.5,[0.5 1.5],'--','color',line_colors(2,:),'linewidth',1.5)
    axis off
    colorbar(ax(1),'NorthOutside')
    for which_pc = 1:3
        ax(which_pc+1) = subplot(4,1,which_pc+1);
        plot(1:length(td_hold),getSig(td_hold,{'M1_pca',which_pc}),'-k')
        hold on
        scatter(1:length(td_hold),getSig(td_hold,{'M1_pca',which_pc}),[],cat(1,td_hold.lambda),'filled')
        ylabel(sprintf('Hold activity in PC %d',which_pc))
        set(gca,'box','off','tickdir','out')
        colormap(viridis)
    end
    xlabel('Trial number')
    linkaxes(ax,'x')
    sgtitle(sprintf('%s %s M1 Hold Activity with \\lambda changes and failures',file_query.monkey,file_query.date))
    
    % Try to fit a linear model between some dimension of hold activity and the previous lambda
    td_hold(1).prev_lambda = td_hold(1).lambda;
    for trialnum = 2:length(td_hold)
        td_hold(trialnum).prev_lambda = td_hold(trialnum-1).lambda;
    end
    [td_hold,lambda_model] = getModel(td_hold,struct(...
        'model_type','linmodel',...
        'model_name','lambda_pred',...
        'in_signals',{{'M1_pca',1:16}},...
        'out_signals',{{'prev_lambda',1}}));
    
    point_color = repmat([0 0 0],length(td_hold),1);
    point_color(lambda_changes+2,:) = repmat([0 1 0],length(lambda_changes),1); % hold activity would change only after 
    figure('defaultaxesfontsize',10)
    scatter(cat(1,td_hold.prev_lambda),cat(1,td_hold.linmodel_lambda_pred),[],point_color,'filled')
    hold on
    plot([1.5 4.5],[1.5 4.5],'--k')
    xlabel('Actual previous \lambda')
    ylabel('Predicted previous \lambda (from hold activity)')
    set(gca,'box','off','tickdir','out')
    axis equal
    title(sprintf('%s %s M1 hold activity previous \\lambda prediction',file_query.monkey,file_query.date))
    
    % Fit a model between hold activity and number of failures in the last K trials
    trial_history = 3;
    is_failure = (cat(1,td_hold.result)=='F');
    num_prev_failures = conv(is_failure,ones(trial_history,1),'full');
    num_prev_failures = num2cell(num_prev_failures(1:length(is_failure)));
    [td_hold.num_prev_failures] = deal(num_prev_failures{:});
    [td_hold,failure_model] = getModel(td_hold,struct(...
        'model_type','linmodel',...
        'model_name','prev_failures',...
        'in_signals',{{'M1_pca',1:16}},...
        'out_signals',{{'num_prev_failures',1}}));
    
    point_color = repmat([0 0 0],length(td_hold),1);
    figure('defaultaxesfontsize',10)
    scatter(cat(1,td_hold.num_prev_failures),cat(1,td_hold.linmodel_prev_failures),[],'k','filled')
    hold on
    plot([0 trial_history],[0 trial_history],'--k')
    xlabel('Actual number of previous failures')
    ylabel('Predicted previous failures (from hold activity)')
    set(gca,'box','off','tickdir','out')
    axis equal
    title(sprintf('%s %s M1 hold activity previous failures prediction',file_query.monkey,file_query.date))
    
%% look for behavioral subspace
%     [~,td_co] = getTDidx(td,'task','CO');
%     [~,td_cst] = getTDidx(td,'task','CST');
%     
%     td_cst = trimTD(td_cst,struct(...
%         'idx_start',{{'idx_cstStartTime',100}},...
%         'idx_end',{{'idx_cstEndTime',0}},...
%         'remove_short',true));
%     
%     td_co = trimTD(td_co,struct(...
%         'idx_start',{{'idx_goCueTime',100}},...
%         'idx_end',{{'idx_goCueTime',500}},...
%         'remove_short',true));
%     
%     [td_cst,cst_pca_info] = dimReduce(td_cst,struct('algorithm','pca','signals','M1_spikes'));
%     [td_co,co_pca_info] = dimReduce(td_co,struct('algorithm','pca','signals','M1_spikes'));
%     
%     td_cst = binTD(td_cst,50);
%     td_co = binTD(td_co,50);
%     
%     [td_cst,cst_hand_vel_decoder] = getModel(td_cst,struct(...
%         'model_type','linmodel',...
%         'model_name','cst_hand_vel_decoder',...
%         'in_signals',{{'M1_pca',1:16}},...
%         'out_signals',{{'hand_vel',1:2}}));
%     
%     [td_co,co_hand_vel_decoder] = getModel(td_co,struct(...
%         'model_type','linmodel',...
%         'model_name','co_hand_vel_decoder',...
%         'in_signals',{{'M1_pca',1:16}},...
%         'out_signals',{{'hand_vel',1:2}}));
%     
%     td_cst = getModel(td_cst,co_hand_vel_decoder);
%     td_co = getModel(td_co,cst_hand_vel_decoder);
%     
%     task_to_plot={'co','cst'};
%     figure
%     for train_task = 1:length(task_to_plot)
%         for test_task = 1:length(task_to_plot)
%             subplot(2,2,(train_task-1)*length(task_to_plot)+test_task)
%             eval(sprintf('test_td=td_%s',task_to_plot{test_task}))
%             
%             scatter(getSig(test_td,{'hand_vel',1}),getSig(test_td,{sprintf('linmodel_%s_hand_vel_decoder',task_to_plot{train_task}),1}),[],'k','filled')
%             xlabel(sprintf('Actual hand velocity in %s',task_to_plot{test_task}))
%             ylabel(sprintf('Predicted hand velocity by %s model',task_to_plot{train_task}))
%             axis equal
%             set(gca,'box','off','tickdir','out','xlim',[-250 350],'ylim',[-250 350])
%         end
%     end
    