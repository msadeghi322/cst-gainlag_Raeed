%%    
    if ispc
        dataroot = 'D:\OneDrive - Northeastern University\Action Lab\01 Projects\Batista Collaboration\00 CST\01 Experiment\data\FordEarl';
    end

    file_info = dir(fullfile(dataroot,'*COCST*'));
    filenames = horzcat({file_info.name})';
    
%% load data
    filenum = 1;
    td_preproc = load_clean_cst_data(fullfile(dataroot,filenames{filenum}));
    
%% compare gpfa to pca
    td = td_preproc;
    
    [~,td] = getTDidx(td,'task','CST','result','R');
    td = trimTD(td,{'idx_cstStartTime',-100},'idx_cstEndTime');

%     % trim to time around go cue
%     start_time = -0.4;
%     end_time = 0.4;
%     td = trimTD(td,struct(...
%         'idx_start',{{'idx_goCueTime',floor(start_time/td(1).bin_size)}},...
%         'idx_end',{{'idx_goCueTime',floor(end_time/td(1).bin_size)}},...
%         'remove_short',true));
    
    num_dims = 10;
    kernSD = 0.03;
    bin_w = 0.02;
    
    fprintf('Running smooth pca...\n')
    td_binned = binTD(td,bin_w/td(1).bin_size);
    td_binned = smoothSignals(td_binned,struct('signals','M1_spikes','width',kernSD,'calc_rate',false,'field_extra','_smooth'));
    td_binned = dimReduce(td_binned,struct(...
        'algorithm','pca',...
        'signals','M1_spikes_smooth',...
        'num_dims',num_dims));
    
    fprintf('Running GPFA...\n')
    td = runGPFA(td,struct(...
        'arrays','M1',...
        'method','gpfa',...
        'xDim',num_dims,...
        'kernSD',kernSD,...
        'bin_w',bin_w));

    [td_binned.M1_gpfa] = deal(td.M1_gpfa);
    
%% Plot GPFA and PCA

dim_red_models = {'smooth_pca','gpfa'};
trial_id = 159;
trial_to_plot = getTDidx(td_binned,'task','CST','trial_id',trial_id);

fig = figure();
for modelnum = 1:length(dim_red_models)
    for compnum = 1:num_dims
        model_trace = td_binned(trial_to_plot).(sprintf('M1_%s',dim_red_models{modelnum}));
        
        subplot(num_dims,length(dim_red_models),(compnum-1)*length(dim_red_models)+modelnum)
        plot(td_binned(trial_to_plot).bin_size*(1:length(model_trace)),model_trace(:,compnum),'-k')
        set(gca,'tickdir','out','box','off','xlim',[0,6],'xtick',[],'ylim',[-1,1])
        
        if modelnum==1
            ylabel(sprintf('Comp %d',compnum))
        end
        
        if compnum==1
            title(dim_red_models{modelnum})
        end
        
        if compnum==num_dims
            xlabel('Time into trial (s)')
            set(gca,'xtick',0:6)
        end
        
    end
end

%% plot gpfa traces
figure
plot_traces(td_binned,struct(...
    'signals',{{'M1_gpfa',1:3}},...
    'trials_to_plot',getTDidx(td_binned,'task','CST','trial_id',159)))

%% compare CO to CST on PCA and GPFA

lambda_to_use = 3.3;

figure
for modelnum = 1:length(dim_red_models)
    % split data
    [~,td_co] = getTDidx(td_binned,'task','CO');
    [~,td_cst] = getTDidx(td_binned,'task','CST');
    
%     % apply pca individually to CST and CO
%     td_co = dimReduce(td_co,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
%     td_cst = dimReduce(td_cst,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));
    
    % plot out neural traces
    dirs = unique(cat(1,td_co.tgtDir));
    dir_colors = linspecer(length(dirs));
    subplot(1,length(dim_red_models),modelnum)
    for dirnum = 1:length(dirs)
        trial_idx = getTDidx(td_co,'task','CO','tgtDir',dirs(dirnum));
        plot_traces(td_co,struct(...
            'signals',{{sprintf('M1_%s',dim_red_models{modelnum}),1:3}},...
            'trials_to_use',trial_idx,...
            'trials_to_plot',trial_idx(randperm(length(trial_idx),10)),...
            'plot_dim',3,...
            'color',dir_colors(dirnum,:)))
        hold on
    end
    trial_idx = getTDidx(td_cst,'task','CST','lambda',lambda_to_use,'trial_id',159);
    plot_traces(td_cst,struct(...
        'signals',{{sprintf('M1_%s',dim_red_models{modelnum}),1:3}},...
        'trials_to_use',trial_idx,...
        'trials_to_plot',trial_idx(randperm(length(trial_idx),1)),...
        'plot_dim',3,...
        'color',[0 0 0]))
    
    title(dim_red_models{modelnum})
    
    set(gca,'box','off','tickdir','out')
    axis equal
    axis off
end