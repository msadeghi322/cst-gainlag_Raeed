% look at CST neural trajectories

%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% Loop through files
lambda_to_use = 3.3;
tangling_timestep = 0.010;
num_dims = 8;
filetic = tic;
for filenum = 27%length(filenames)
    td = load_clean_cst_data(fullfile(dataroot,'library',filenames{filenum}));
    
    % Make sure we have CST trials
    if isempty(td)
        fprintf('Incomplete dataset for file %d\n',filenum)
        continue
    end
    
    if isempty(td(1).M1_unit_guide)
        fprintf('Skipping file %d because no spike data...\n',filenum)
        continue
    end
    
    % get rid of unsorted neurons
    bad_units = td(1).M1_unit_guide(:,2)<=1;
    % for file 31 only! remove neurons 8,2 and 64,2
    if contains(td(1).date_time,'2018/06/26')
        corr_units = [8 2;64 2];
        bad_units = bad_units | ismember(td(1).M1_unit_guide,corr_units,'rows');
    end
    for trialnum = 1:length(td)
        td(trialnum).M1_spikes = td(trialnum).M1_spikes(:,~bad_units);
        td(trialnum).M1_unit_guide = td(trialnum).M1_unit_guide(~bad_units,:);
    end
    
    td = removeBadNeurons(td,struct(...
        'do_fr_check',true,...
        'min_fr',0.1,...
        'fr_window',{{'idx_cstStartTime',0;'idx_cstEndTime',0}},...
        'calc_fr',true,...
        'use_trials',getTDidx(td,'lambda',lambda_to_use)));
    
    % smooth data
    td = smoothSignals(td,struct('signals','M1_spikes','width',0.075,'calc_rate',true));
    td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    td = dimReduce(td,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
    td = getDifferential(td,struct('signals','M1_pca','alias','M1_pca_diff'));
    
    % split data
    [~,td_co] = getTDidx(td,'task','CO');
    [~,td_cst] = getTDidx(td,'task','CST');
    
    % trim TD to only CST portion
    td_co = trimTD(td_co,'idx_goCueTime','idx_otHoldTime');
    td_cst = trimTD(td_cst,'idx_cstStartTime','idx_cstEndTime');
    
    % bin at 10 ms (or whatever timestep we choose)
    td_co = binTD(td_co,floor(tangling_timestep/td_co(1).bin_size));
    td_cst = binTD(td_cst,floor(tangling_timestep/td_cst(1).bin_size));
    
    % get baseline tangling from CO
    q_co = calc_tangling(cat(1,td_co.M1_pca),cat(1,td_co.M1_pca_diff));
    start_idx = 1;
    for trialnum = 1:length(td_co)
        end_idx = size(td_co(trialnum).M1_pca,1)+start_idx-1;
        td_co(trialnum).M1_tangling = q_co(start_idx:end_idx);
        start_idx = end_idx+1;
    end
    
    % check tangling over all trials in a particular lambda
    [~,td_lambda] = getTDidx(td_cst,'lambda',lambda_to_use);
    q_cst = calc_tangling(cat(1,td_lambda.M1_pca),cat(1,td_lambda.M1_pca_diff));
    start_idx = 1;
    for trialnum = 1:length(td_lambda)
        end_idx = size(td_lambda(trialnum).M1_pca,1)+start_idx-1;
        td_lambda(trialnum).M1_tangling = q_cst(start_idx:end_idx);
        start_idx = end_idx+1;
    end
    
    % find limits of tangling
    tangling_clims = prctile(log10(cat(1,q_co,q_cst)),[1 90]);
    
    % make plots to show tangling
    co_plot_params = [...
        struct(...
            'plot_loc',1,...
            'plot_fcn',@(trial) scatter(...
                trial.hand_pos(:,1),...
                -trial.hand_pos(:,2)+470,...
                [],log10(trial.M1_tangling),'filled'),...
            'xlabel','',...
            'ylabel','',...
            'title','Hand movement',...
            'gca_params',struct(...
                'box','off',...
                'tickdir','out',...
                'xlim',[-60 60],...
                'ylim',[-120 120],...
                'clim',tangling_clims,...
                'DataAspectRatio',[1 1 1],...
                'DataAspectRatioMode','manual')),...
        struct(...
            'plot_loc',2,...
            'plot_fcn',@(trial) scatter3(...
                trial.M1_pca(:,1),...
                trial.M1_pca(:,2),...
                trial.M1_pca(:,3),...
                [],log10(trial.M1_tangling),'filled'),...
            'xlabel','',...
            'ylabel','',...
            'title','M1 PCA',...
            'gca_params',struct(...
                'box','off',...
                'tickdir','out',...
                'clim',tangling_clims,...
                'DataAspectRatio',[1 1 1],...
                'DataAspectRatioMode','manual'))];
    handle = make_interactive_td_plot(td_co,...
        struct(...
            'shape',[1 2],...
            'colormap',viridis,...
            'suptitle_fcn',...
                @(trial) sprintf('%s %s Trial ID: %d',...
                trial.monkey,...
                trial.date_time,...
                trial.trial_id),...
            'gcf_params',struct(...
                'defaultaxesfontsize',18,...
                'position',[1544 397 813 450])),...
        co_plot_params);
    
%     h = plot_interactive_cst_phase(td_lambda,struct(...
%         'cursor_sig',{{'cursor_pos',1}},...
%         'hand_sig',{{'hand_pos',1}},...
%         'color_sig',{{'M1_tangling',1}}));

    fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
end