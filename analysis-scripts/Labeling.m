
%% Labeling trials based on control strategy

clear
clc
close all
set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultAxesFontSize',14)
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
set(groot,'defaultAxesBox','off')
set(0, 'DefaultFigureRenderer', 'painters');


%% ************************* Important settings ***************************

% --------------  Choose subject --------------
% 1: Ford, 2: Earl

CriticalLambda = [3.5 , 2];

% -------------- Classification threshold -----
% This is the ratio of cursor RMS_pos over RMS_vel. Larger ratio means
% stronger velocity control strategy. By increasing the threshold
% the classification becomes more strict in identifying a trial as velocity
% control. This threshold is determined based on the simulations and should
% have a minimum of 0.5 (reasonable range [.7--1]). 
Threshold_low = .5;
Threshold_up = .85;


% Would you like to visualize individual trials with the labels for each
% lambda value?
IndividualPlots = 1; % 1:Yes  ,   0:No



% The labels are stored in the following variables. The trial numbers are
% counted after concatination of original data files.

%   PosContTrlIdx: 
%   contains the trial numbers for trials with position control

%   VelContTrlIdx: 
%   contains the trial numbers for trials with velocity control


% ________________________________________________________________________

FilePath = '/Users/Shared/Previously Relocated Items/Security/00 Postdoc at Northeastern/01 Projects/Batista Collaboration/00 CST /cst-gainlag_Raeed/All data/';
FileName = {'Ford_20180618_COCST_TD'
            'Ford_20180626_COCST_TD'
            'Ford_20180627_COCST_TD'
            'Earl_20190716_COCST_TD'
            'Earl_20190702_COCST_TD'};


% if SubjectID==1
%     FileName = {'data\FordEarl\All data\Ford_20180618_COCST_TD'
%                 'data\FordEarl\All data\Ford_20180626_COCST_TD'
%                 'data\FordEarl\All data\Ford_20180627_COCST_TD'};
% elseif SubjectID==2
%     FileName = {'data\FordEarl\All data\Earl_20190716_COCST_TD'
%                 'data\FordEarl\All data\Earl_20190702_COCST_TD'};
% else
%     error('Subject ID is incorrect')
% end
        
        

%% ************************* Concatinate Data *****************************
tic
Variables = {'monkey'
             'task'
             'date'
             'trial_id'
             'lambda'
             'result'
             'idx_startTime'
             'idx_goCueTime'
             'idx_endTime'
             'hand_pos'
             'cursor_pos'
             } ;

DataStruct = struct;
DataTable = table;
trcnt = 0; % counting trials
for i=1:length(FileName)
    tmp = sprintf('%s%s',FilePath,FileName{i});
    D = load(tmp);
    D = D.trial_data;
    T = struct2table(D);
    dt = T.date_time{1}(1:10);  % only extract the date
    T.date = cell(size(T,1),1); % add date as new variable
    T.date(:) = {dt};
    ind = ismember(  T.Properties.VariableNames' , Variables  );
    DataTable = [DataTable ; T(:,ind)];

end

LoadTime = toc ;

%% Calculate RMS
TT = DataTable;

% Data filtering
fc  = 6;      % cut-off frequency
fps = 1000;
[B,A] = butter(6,2*fc/fps,'low');
delta = .001;
% Initialization
lngth = size(TT,1);
TT.HandRMS_p     = NaN(lngth,1);
TT.CursRMS_p     = NaN(lngth,1);
TT.HandRMS_v     = NaN(lngth,1);
TT.CursRMS_v     = NaN(lngth,1);
TT.CursMeanPos   = NaN(lngth,1);
TT.CursMeanVel   = NaN(lngth,1);
TT.RMSRatio      = NaN(lngth,1);
TT.ControlPolicy = cell(lngth,1); TT.ControlPolicy(:)={'Null'};

tic
for i = 1:lngth
    
    % if trial is center out or failed, ignore
    if contains(TT.result(i),'F') || contains(TT.task(i),'CO')
        continue
    end
    
    p1   = TT.idx_goCueTime(i)+50;
    p2   = TT.idx_endTime(i)-500;
    C_p  = filtfilt(B,A, TT.cursor_pos{i}(p1:p2,1) )/1000;
    H_p  = filtfilt(B,A, TT.hand_pos{i}(p1:p2,1) )/1000;
    % centre the cursor and hand position based on cursor first sample
    ref =  C_p(1);
    H_p = H_p - ref; % cursor position is the reference
    C_p = C_p - ref;
    Time = (1:length(C_p))'*delta;
    
    ii   = isnan(H_p);
    H_p(ii)=[];
    C_p(ii)=[];
    Time(ii)=[];
    H_v  = diff(H_p)/delta; H_v=filtfilt(B,A,[H_v;H_v(end)]);
    C_v  = diff(C_p)/delta; C_v=filtfilt(B,A,[C_v;C_v(end)]);
    
   
    TT.HandRMS_p(i,1)   = rms(H_p);
    TT.CursRMS_p(i,1)   = rms(C_p);
    TT.HandRMS_v(i,1)   = rms(H_v);
    TT.CursRMS_v(i,1)   = rms(C_v);
    TT.CursMeanPos(i,1) = mean(C_p);
    TT.CursMeanVel(i,1) = mean(C_v);
    

    
end
lapsedtime = toc;

%% Classification

% parameters: 
% ConfidenceLevel: The level of confidence in classification
% LambdaThresh: Upper limit for lambda value to be considered in the
% analysis

ConfidenceLevel = 55;
LambdaThresh = 5;% max( CriticalLambda );

SVMModel = load('SVMModel');
SVMModel = SVMModel.SVMModel;
CompactSVM = compact(SVMModel);
CompactSVM = fitPosterior(CompactSVM, SVMModel.X, SVMModel.Y);

DataTable = TT;
idx = contains(TT.result,'R') & contains(TT.task,'CST') & TT.lambda<=LambdaThresh;
y1 = TT.CursRMS_p(idx)*100;
y2 = TT.CursRMS_v(idx)*100;
[PredictedLabel,Score] = predict(CompactSVM,[y1 , y2]);
ii = PredictedLabel==0 & Score(:,1)>ConfidenceLevel/100;
jj = PredictedLabel==1 & Score(:,2)>ConfidenceLevel/100;
kk = (~ii) & (~jj);

idxN = find(idx==1);
TT.ControlPolicy(idxN(ii)) = {'Position'};
TT.ControlPolicy(idxN(jj)) = {'Velocity'};
TT.ControlPolicy(idxN(kk)) = {'Hybrid'};



%% Save the table 
Labels = TT;
Labels.idx_startTime = [];
Labels.idx_goCueTime = [];
Labels.idx_endTime = [];
Labels.hand_pos = [];
Labels.cursor_pos = [];

save(sprintf('Labels_conf_%d',ConfidenceLevel),'Labels');




%% Classification
GroupColor = [.2,.2,.2
              .5,.4,0
              0,.6,.6 ];
Color2 = [.2 .2 .2]; 
Color3 = [.2,.2,0 ; 0,.2,.2];
FontSize = 11;
MarkerSize = 12;
LineWidth = 1.5;
Mnk = {'Ford','Earl'};


SVMModel   = load('SVMModel');
SVMModel   = SVMModel.SVMModel;
CompactSVM = compact(SVMModel);
CompactSVM = fitPosterior(CompactSVM, SVMModel.X, SVMModel.Y);
Thresh = 3;

ffg = figure(3);
clf;
sb1 = 3;
sb2 = 6; %length(SList);
for ID=1:length(Mnk)
    id = contains(Labels.monkey , Mnk{ID}) & ~isnan(Labels.lambda);
    D = Labels;
    Lmbda = D.lambda;
    ll = id==1 & Lmbda <= Thresh*CriticalLambda(ID);
    
    y1 = Labels.CursRMS_p(ll)*100;
    y2 = Labels.CursRMS_v(ll)*100;
    
    subplot(sb1,sb2,ID)
    hold all
    
    % classified indexes
    [PredictedLabel,Score] = predict(CompactSVM,[y1 , y2]);
    ii = PredictedLabel==0 & Score(:,1)>ConfidenceLevel/100;
    jj = PredictedLabel==1 & Score(:,2)>ConfidenceLevel/100;
    kk = (~ii) & (~jj);
    
    % some example trials
    i1 = PredictedLabel==0 & Score(:,1)>.8; % high conf pos trials
    i2 = PredictedLabel==1 & Score(:,2)>.8; % high conf vel trials
    pos_exp_x = y1(i1);  [pos_exp_x,j1] = sort(pos_exp_x);
    pos_exp_y = y2(i1);  pos_exp_y = pos_exp_y(j1);
    vel_exp_x = y1(i2);  [vel_exp_x,j2] = sort(vel_exp_x);
    vel_exp_y = y2(i2);  vel_exp_y = vel_exp_y(j2);
    
    h1 = plot(y1(kk),y2(kk),'o','MarkerEdgeColor',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'markersize',MarkerSize-9,'Linewidth',1);
    h2 = plot(y1(ii),y2(ii),'o','MarkerEdgeColor',GroupColor(2,:),'MarkerFaceColor',[1 1 1],'markersize',MarkerSize-9,'Linewidth',1);
    h3 = plot(y1(jj),y2(jj),'o','MarkerEdgeColor',GroupColor(3,:),'MarkerFaceColor',[1 1 1],'markersize',MarkerSize-9,'Linewidth',1);
    plot(pos_exp_x(end),pos_exp_y(end),'s','MarkerEdgeColor',Color3(1,:),'MarkerFaceColor',[1 1 1],'markersize',MarkerSize-3,'Linewidth',1);
    plot(vel_exp_x(end),vel_exp_y(end),'^','MarkerEdgeColor',Color3(2,:),'MarkerFaceColor',[1 1 1],'markersize',MarkerSize-3,'Linewidth',1);
    plot(pos_exp_x(end),pos_exp_y(end),'.','MarkerEdgeColor',GroupColor(2,:),'MarkerFaceColor',[1 1 1],'markersize',MarkerSize-6,'Linewidth',1);
    plot(vel_exp_x(end),vel_exp_y(end),'.','MarkerEdgeColor',GroupColor(3,:),'MarkerFaceColor',[1 1 1],'markersize',MarkerSize-6,'Linewidth',1);
    PosProb(ID) = mean( Score(~kk,1) );
    VelProb(ID) = mean( Score(~kk,2) );
    text(.5,4,sprintf('P(pos.)= %.2f',PosProb(ID)))
    xlim([0,5])
    ylim([0,5])
    set(gca,'fontsize',FontSize)
    if ID==1
        ylabel(sprintf('RMS Vel (cm/s)'))
        legend('Hybrid','Pos','Vel')
    end
    xlabel(sprintf('RMS Pos (cm)'))
    title(Mnk{ID})
    
    
    
end


subplot(sb1,sb2,4)
hold all
y = PosProb;
x = 3;
plot(x,y,'p','color',GroupColor(1,:),'MarkerFaceColor',GroupColor(1,:),'markersize',MarkerSize-5);
xlim([.5 3.5]);
ylim([0,1])
plot([0,4],[.5,.5],':k');
set(gca,'xtick',3,'xticklabel',{'Monkey'})
%xlabel('Group')
ylabel('P(position control)')
set(gca,'fontsize',FontSize)







return
%%

% Extract some useful information
D = DataStruct;
delta = .001;
TrialN = length(D);
for tr=1:TrialN
    Lambda(tr,1) = D(tr).lambda;
    Result(tr,1) = D(tr).result;
    StartTimeIdx(tr,1) = D(tr).idx_startTime(end);
    GoCueTimeIdx(tr,1) = D(tr).idx_goCueTime(end);
    TaskType{tr,1} = D(tr).task;
end

trl_cst = find(~isnan(Lambda));
ii = isnan(Lambda);
LambdaList = unique(Lambda(~ii));
LambdaCount = zeros(size(LambdaList));
for lm=1:length(LambdaList)
    jj = Lambda==LambdaList(lm);
    LambdaCount(lm,1) = sum(jj);
    RW = Result(jj);
    i_succ = RW =='R';
    SuccessRate(lm,1) = sum(i_succ)/sum(jj); 
end




%% *********************** Calculate the metrics **************************

% Data filtering
fc  = 6;      % cut-off frequency
fps = 1000;
[B,A] = butter(6,2*fc/fps,'low');

% Initialization
LN = length(LambdaList);
TN = max(LambdaCount);
HandRMS_p     = NaN(TN,LN);
CursRMS_p     = NaN(TN,LN);
HandRMS_v     = NaN(TN,LN);
CursRMS_v     = NaN(TN,LN);
CursMeanPos   = NaN(TN,LN);
CursMeanVel   = NaN(TN,LN);
TrialNumber   = NaN(TN,LN);
ControlStrategy = NaN()

for i = 1:length(LambdaList)
    idj = find(Lambda==LambdaList(i));
    cc=0;
    for j=1:length(idj)
        tr = idj(j);
        if contains(D(tr).result,'F')
            continue
        end
        
        p1   = D(tr).idx_goCueTime+50;
        p2   = D(tr).idx_endTime-500;
        C_p  = filtfilt(B,A, D(tr).cursor_pos(p1:p2,1) )/1000;
        H_p  = filtfilt(B,A,D(tr).hand_pos(p1:p2,1) )/1000;
        % centre the cursor and hand position based on cursor first sample
        ref =  C_p(1);
        H_p = H_p - ref; % cursor position is the reference
        C_p = C_p - ref;
        
        Time = (1:length(C_p))'*delta;
        ii   = isnan(H_p);
        H_p(ii)=[];
        C_p(ii)=[];
        Time(ii)=[];
        H_v  = diff(H_p)/delta; H_v=filtfilt(B,A,[H_v;H_v(end)]);
        C_v  = diff(C_p)/delta; C_v=filtfilt(B,A,[C_v;C_v(end)]);
        
        TrialNumber(j,i) = tr;
        HandRMS_p(j,i) = rms(H_p);
        CursRMS_p(j,i) = rms(C_p);
        HandRMS_v(j,i) = rms(H_v);
        CursRMS_v(j,i) = rms(C_v);
        CursMeanPos(j,i) = mean(C_p);
        CursMeanVel(j,i) = mean(C_v);
    end
end




%% Classification based on RMS ratio: RMS_pos / RMS_vel


figure(300)
sb1 = 4;
sb2 = 5;
clf

VelContTrlIdx_mtrx = NaN(size(TrialNumber));
PosContTrlIdx_mtrx = NaN(size(TrialNumber));
for i=1:size(CursRMS_p,2)
   
    X = CursRMS_p(:,i); 
    Y = CursRMS_v(:,i); 
    Z = X./Y;
    jj = Z>Threshold_low;
    VelContTrlIdx_mtrx(jj , i) = TrialNumber(jj,i);
    PosContTrlIdx_mtrx(~jj ,i) = TrialNumber(~jj,i);
    
    ii=~isnan(X);
    X = X(ii);
    Y = Y(ii);
    ZZ = X./Y;
    
    subplot(sb1,sb2,i)
    hold all
    plot(ZZ,'.k');
    plot(0*ZZ+Threshold_low,'--r')
    ylim([0 2]);
    xlabel('trial')
    ylabel('RMS Ratio')
    set(gca,'fontsize',11)
    text(1,1.9,sprintf('\\lambda = %.2f',LambdaList(i)),'fontsize',12)
    
end

VelContTrlIdx = sort( VelContTrlIdx_mtrx(:) ); VelContTrlIdx(isnan(VelContTrlIdx)) = [];
PosContTrlIdx = sort( PosContTrlIdx_mtrx(:) ); PosContTrlIdx(isnan(PosContTrlIdx)) = [];


%%
CL = [.6,.5,0
      0,.7,.6];
  
figure(400)
sb1 = 3;
sb2 = 4;
clf

IND = {[1:size(CursMeanPos,2)];
       [1:ceil(size(CursMeanPos,2)/2)];
       [ceil(size(CursMeanPos,2)/2):size(CursMeanPos,2)]};
PData = {CursMeanPos , CursMeanVel , CursRMS_p , CursRMS_v };


for lv = 1:length(IND)
    
    ind = IND{lv};
    
    %--------- Mean ----------
    
    X = CursMeanPos(:,ind); X=X(:);
    Y = CursMeanVel(:,ind); Y=Y(:);
    Z_p = PosContTrlIdx_mtrx(:,ind); Z_p = Z_p(:); % trial idx for position control
    Z_v = VelContTrlIdx_mtrx(:,ind); Z_v = Z_v(:); % trial idx for position control
    
    iip=~isnan(Z_p); % find pos control trials
    Xp = X(iip);
    Yp = Y(iip);
    
    iiv=~isnan(Z_v); % find vel control trials
    Xv = X(iiv);
    Yv = Y(iiv);
    
    ii=~isnan(X); % all trials that are not NaN
    XX = X(ii);
    YY = Y(ii);
    [rho,~] = corr(XX,YY);
    
    subplot(sb1,sb2,lv)
    hold all
    hp = plot(Xp,Yp,'.','color',CL(1,:));
    hv = plot(Xv,Yv,'.','color',CL(2,:));
    plot([-.05 .05],[0,0],':k');
    plot([0,0],[-.01 .01],':k');
    ylim([-.01 .01]);
    xlim([-.05 .05])
    xlabel('Mean pos.')
    ylabel('Mean vel.')
    set(gca,'fontsize',11)
    text(-.045,.0075,sprintf('R = %.2f',rho),'fontsize',12)
    title('All trials')
end
legend([hp,hv],'PosCont','VelCont','location','best')



for lv = 1:length(IND)
    
    ind = IND{lv};
    
    %--------- Mean ----------
    subplot(sb1,sb2,lv+sb2)
    hold all
    X = CursRMS_p(:,ind); X=X(:);
    Y = CursRMS_v(:,ind); Y=Y(:);
    Z_p = PosContTrlIdx_mtrx(:,ind); Z_p = Z_p(:); % trial idx for position control
    Z_v = VelContTrlIdx_mtrx(:,ind); Z_v = Z_v(:); % trial idx for position control
    
    iip=~isnan(Z_p); % find pos control trials
    Xp = X(iip);
    Yp = Y(iip);
    
    iiv=~isnan(Z_v); % find vel control trials
    Xv = X(iiv);
    Yv = Y(iiv);
    
    ii=~isnan(X); % all trials that are not NaN
    XX = X(ii);
    YY = Y(ii);
    
    hp = plot(Xp,Yp,'.','color',CL(1,:));
    hv = plot(Xv,Yv,'.','color',CL(2,:));
    plot([-.05 .05],[0,0],':k');
    plot([0,0],[-.01 .01],':k');
    ylim([0 .06]);
    xlim([0 .04])
    xlabel('RMS pos.')
    ylabel('RMS vel.')
    set(gca,'fontsize',11)
    
end
legend([hp,hv],'PosCont','VelCont','location','best')



%%

if ~IndividualPlots
return
end

for i = 1:length(LambdaList)
    
    figure(i)
    sb1 = 7;
    sb2 = 9;
    clf
    
    
    idj = find(Lambda==LambdaList(i));
    cc=0;
    for j=1:length(idj)
        
        tr = idj(j);
        if contains(D(tr).result,'F')
            continue
        end
        
        Flag_v = ismember(tr,VelContTrlIdx); % is this trial vel control?
        
        p1   = D(tr).idx_goCueTime+50;
        p2   = D(tr).idx_endTime-500;
        H_p  = filtfilt(B,A,D(tr).hand_pos(p1:p2,1) )/1000;
        C_p  = filtfilt(B,A, D(tr).cursor_pos(p1:p2,1) )/1000;
        % centre the cursor and hand position based on cursor first sample
        ref =  C_p(1);
        H_p = H_p - ref; % cursor position is the reference
        C_p = C_p - ref;
        
        Time = (1:length(C_p))'*delta;
        ii   = isnan(H_p);
        H_p(ii)=[];
        C_p(ii)=[];
        Time(ii)=[];
        
        
        H_v  = diff(H_p)/delta; H_v=filtfilt(B,A,[H_v;H_v(end)]);
        C_v  = diff(C_p)/delta; C_v=filtfilt(B,A,[C_v;C_v(end)]);
        L    = D(tr).lambda;
        
        cc=cc+1;
        if cc>sb1*sb2; continue; end
        subplot(sb1,sb2,cc)
        hold all
        H1 = plot(Time,C_p,'b','linewidth',1.5);
        H2 = plot(Time,H_p,'r','linewidth',1.5);
        plot(Time(end)*[0,1],[0,0],':k')
        if Flag_v
            plot([0,6],[.05 .05],'-','color',CL(2,:),'linewidth',3);
        else
            plot([0,6],[.05 .05],'-','color',CL(1,:),'linewidth',3);
        end
        
        
        if cc==1
            legend([H1,H2],'Cursor','Hand','location','best')
            title(sprintf('\\lambda = %.3f',L));
            xlabel('Time (s)')
            ylabel('Position')
        end
        ylim([-.05 .05])
        xlim(6*[0,1])
        set(gca,'FontSize',8)
        
        
    end
    
    
    
    
end



