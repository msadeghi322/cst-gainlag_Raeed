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
    MonkeyNames = {'Earl','Ford'};
    

    %% Data extraction:
    
    % To make life easier, only extract data from CST trials that are
    % successful. Eliminate all the others. 
    
    if ~exist(sprintf('%s/DataTableAll.mat',dataroot),'file')
        DataTable = ms_cst_datatable_all(dataroot, filenames);
        save(sprintf('%s/DataTableAll',dataroot), 'DataTable' , '-v7.3');
    else
        D = load(sprintf('%s/DataTableAll',dataroot));
        DataTable = D.DataTable;
    end
    
    
    %% Pick monkey and date
    
    MonkeyID = 2;           % 1:Earl 2:Ford
    Session = 3;            % 1,2: Earl , 3-5:Ford 
    num_dims = 8;           % number of dimensions for pca
    softnorm = true;        %  
    smoothsigs = true;
    
    
    
    idx = contains(DataTable.monkey, MonkeyNames(MonkeyID)) & ...
        contains( DataTable.date , fileDates(Session) );
    td_preproc = table2struct(DataTable(idx,:));
    td = td_preproc;
    
    
    
    %% Smoothing the data
    
    
    
    if smoothsigs
        td = smoothSignals(td,struct('signals','M1_spikes','width',0.050,'calc_rate',true));
    end

    % trim TD only for successful trials
    startSample = 500;
    endSample   = 5500;
    td = trimTD(td,{'idx_goCueTime',startSample},{'idx_goCueTime',endSample});

    % smooth data
    if softnorm
        td = softNormalize(td,struct('signals','M1_spikes','alpha',5));
    end

    
    %%  Dimensionality reduction
    
    % get average for each bin
    td_binned = binTD(td,'average');
    td_binned = dimReduce(td_binned,struct('algorithm','pca','signals','M1_spikes','num_dims',num_dims));

    % split data
    %[~,td_co]  = getTDidx(td_binned,'task','CO');
    %[~,td_cst] = getTDidx(td_binned,'task','CST');
    %[~,td_co] = getTDidx(td_binned,'task','CST','lambda',4.1);
    td_cst = td_binned;
    
    CstTrlId = [td_cst.trial_id]';
    
    
    %% Find trial index for control strategy
    
    ConfidenceLevel = 51; % confidence percentage for classification

    Labels = load(sprintf('Labels_conf_%d',ConfidenceLevel)); Labels = Labels.Labels;
    
    iip = contains(Labels.monkey,MonkeyNames{MonkeyID}) & ...
        contains(Labels.date, fileDates(Session)) & ...
         contains(Labels.ControlPolicy,'Position');
    PosTrlId =  Labels.trial_id(iip);
    
    
    iiv = contains(Labels.monkey,MonkeyNames{MonkeyID}) & ...
        contains(Labels.date, fileDates(Session)) & ...
         contains(Labels.ControlPolicy,'Velocity');
    VelTrlId =  Labels.trial_id(iiv);
    
        
    iih = contains(Labels.monkey,MonkeyNames{MonkeyID}) & ...
        contains(Labels.date, fileDates(Session)) & ...
         contains(Labels.ControlPolicy,'Hybrid');
    HybTrlId =  Labels.trial_id(iih);
    

    
   
    %--------------------- CMD Scaling / PCA
    
    M1 = cat(1,td_cst.M1_spikes);
    Y = pdist(M1);
    Z = squareform(Y);
    [Yc , ic] = cmdscale(Z,5);

    
    figure(1)
    clf
    sb1=3;
    sb2=6;
    MS = 12;
    GroupCl = [.8,.5,.0 ; .0,.7,.7 ; .7,.7,.7];
    subplot(sb1,sb2,1)
    ip = ismember(CstTrlId , PosTrlId);
    iv = ismember(CstTrlId , VelTrlId);
    ih = ismember(CstTrlId , HybTrlId);
    plot3(Yc(ip,3),Yc(ip,2),Yc(ip,1),'.','color',GroupCl(1,:),'markersize',MS); hold on
    plot3(Yc(iv,3),Yc(iv,2),Yc(iv,1),'.','color',GroupCl(2,:),'markersize',MS);
    plot3(Yc(ih,3),Yc(ih,2),Yc(ih,1),'.','color',GroupCl(3,:),'markersize',MS);
    grid
    xlabel('Config 3')
    ylabel('Config 2')
    zlabel('Config 1')
    title('Multi-Dimensional Scaling')
    
    
    subplot(sb1,sb2,2)
    plot(Yc(ip,2),Yc(ip,1),'.','color',GroupCl(1,:),'markersize',MS); hold on
    plot(Yc(iv,2),Yc(iv,1),'.','color',GroupCl(2,:),'markersize',MS);
    plot(Yc(ih,2),Yc(ih,1),'.','color',GroupCl(3,:),'markersize',MS);
    xlabel('Config 2')
    ylabel('Config 1')
    box off
    axis equal

    subplot(sb1,sb2,3)
    plot(Yc(ip,3),Yc(ip,1),'.','color',GroupCl(1,:),'markersize',MS); hold on
    plot(Yc(iv,3),Yc(iv,1),'.','color',GroupCl(2,:),'markersize',MS);
    plot(Yc(ih,3),Yc(ih,1),'.','color',GroupCl(3,:),'markersize',MS);
    xlabel('Config 3')
    ylabel('Config 1')
    box off
    axis equal

    subplot(sb1,sb2,4)
    plot(Yc(ip,3),Yc(ip,2),'.','color',GroupCl(1,:),'markersize',MS); hold on
    plot(Yc(iv,3),Yc(iv,2),'.','color',GroupCl(2,:),'markersize',MS);
    plot(Yc(ih,3),Yc(ih,2),'.','color',GroupCl(3,:),'markersize',MS);
    xlabel('Config 3')
    ylabel('Config 2')
    legend('Pos','Vel','None')
    box off
    axis equal

    
    
    C = copper(size(Yc,1));
    
    subplot(sb1,sb2,sb2+1)
    S = scatter3(Yc(:,3),Yc(:,2),Yc(:,1),MS*10,C);
    S.Marker = '.';
    xlabel('Config 3')
    ylabel('Config 2')
    zlabel('Config 1')
    
    subplot(sb1,sb2,sb2+2)
    S = scatter(Yc(:,2),Yc(:,1),MS*10,C);
    S.Marker = '.';
    xlabel('Config 2')
    ylabel('Config 1')
    box off
    axis equal

    subplot(sb1,sb2,sb2+3)
    S = scatter(Yc(:,3),Yc(:,1),MS*10,C);
    S.Marker = '.';
    xlabel('Config 3')
    ylabel('Config 1')
    box off
    axis equal

    subplot(sb1,sb2,sb2+4)
    hold all
    S = scatter(Yc(:,3),Yc(:,2),MS*10,C);
    S.Marker = '.';
    xlabel('Config 3')
    ylabel('Config 2')
    box off
    axis equal
    h1 = plot(Yc(1,3),Yc(1,2),'.','color',C(1,:), 'markersize', 10 );
    h2 = plot(Yc(end,3),Yc(end,2),'.','color',C(end,:), 'markersize', 10);
    legend([h1,h2],'Early','Late')
    
    
    %----------------------- T-SNE analysis
    
    M1 = cat(1,td_cst.M1_spikes);
    TrID = cat(1,td_cst.trial_id);
    Y_sne = tsne(M1);
    
    figure(2)
    clf
    sb1=3;
    sb2=6;
    subplot(sb1,sb2,1)
    hold all
    plot(Y_sne(ip,2),Y_sne(ip,1),'.','color',GroupCl(1,:),'markersize',MS)
    plot(Y_sne(iv,2),Y_sne(iv,1),'.','color',GroupCl(2,:),'markersize',MS)
    plot(Y_sne(ih,2),Y_sne(ih,1),'.','color',GroupCl(3,:),'markersize',MS,'linewidth',1.5)
    xlabel('dim2')
    ylabel('dim1')
    axis equal
    title(sprintf('t-sne: [%d  %d]ms', startSample , endSample))
    legend('Pos','Vel','None')
    
    subplot(sb1,sb2,sb2+1)
    hold all
    C = copper(size(Y_sne,1));
    S = scatter(Y_sne(:,2),Y_sne(:,1),10*MS,C);
    S.Marker = '.';
    xlabel('dim2')
    ylabel('dim1')
    axis equal
    h1 = plot(Y_sne(1,2),Y_sne(1,1),'.','color',C(1,:), 'markersize', 10 );
    h2 = plot(Y_sne(end,2),Y_sne(end,1),'.','color',C(end,:), 'markersize', 10);
    legend([h1,h2],'Early','Late')
    
    
    
    
    
    return
    
    
    
    % PCA analysis
    
    figure(3)
    clf
    sb1=2;
    sb2=4;
    
    M1_pca = cat(1,td_cst.M1_pca);
    
    subplot(sb1,sb2,1)
    plot3(M1_pca(ip,3),M1_pca(ip,2),M1_pca(ip,1),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot3(M1_pca(iv,3),M1_pca(iv,2),M1_pca(iv,1),'.','color',GroupCl(2,:),'markersize',16);
    plot3(M1_pca(ih,3),M1_pca(ih,2),M1_pca(ih,1),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('PC3')
    ylabel('PC2')
    zlabel('PC1')
    grid
    axis equal

    subplot(sb1,sb2,2)
    plot(M1_pca(ip,2),M1_pca(ip,1),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot(M1_pca(iv,2),M1_pca(iv,1),'.','color',GroupCl(2,:),'markersize',16);
    plot(M1_pca(ih,2),M1_pca(ih,1),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('PC2')
    ylabel('PC1')
    legend('Pos','Vel','None')
    box off
    axis equal

    
    subplot(sb1,sb2,3)
    plot(M1_pca(ip,3),M1_pca(ip,1),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot(M1_pca(iv,3),M1_pca(iv,1),'.','color',GroupCl(2,:),'markersize',16);
    plot(M1_pca(ih,3),M1_pca(ih,1),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('PC3')
    ylabel('PC1')
    legend('Pos','Vel','None')
    box off
    axis equal

    
    subplot(sb1,sb2,4)
    plot(M1_pca(ip,3),M1_pca(ip,2),'.','color',GroupCl(1,:),'markersize',16); hold on
    plot(M1_pca(iv,3),M1_pca(iv,2),'.','color',GroupCl(2,:),'markersize',16);
    plot(M1_pca(ih,3),M1_pca(ih,2),'.','color',GroupCl(3,:),'markersize',16);
    xlabel('PC3')
    ylabel('PC2')
    legend('Pos','Vel','None')
    box off
    axis equal
    
    
    
    
    
    %% %%%%%%%%%%%%%% Local functions
    
    function DataTable = ms_cst_datatable_all(dataroot, filenames)
    
    % This function finds all the data mat files in a given folder and
    % concatinates them all in a table. The current mat files should be of type
    % structure. The resultant data table only contains the variables that are
    % shared among all the mat files.
    
    DataTable = table;
    for ID=1:length(filenames)
        
        D = load_clean_cst_data(fullfile(dataroot,sprintf('%s',filenames{ID})));
        T = struct2table(D);
        dt = T.date_time{1}(1:10);  % only extract the date
        T.date = cell(size(T,1),1); % add date as new variable
        T.date(:) = {dt};
        idx1 = contains(T.task,'CO')| contains(T.result,'F');
        T(idx1,:) = [];
        ind = ~ismember(DataTable.Properties.VariableNames , T.Properties.VariableNames);
        DataTable(:,ind) = [];
        %ind = ismember(  T.Properties.VariableNames' , Variables  );
        DataTable = [DataTable ; T];
    end
    
    
    end
    
    
    
    
