
function NewTable = RefineLabels (OldTable , Threshold_low , Threshold_up)

% OldTable: the current table

% Threshold_low: range between 0.45 to 0.7: determines the position control
% policy. 

% Threshold_low: range between 0.45 to 0.7: determines the position control
% policy



if nargin<2  % default values
    Threshold_low = .5;
    Threshold_up  = .85;
elseif nargin <3
    Threshold_up  = .85;
end

% Set the threshold
TT = OldTable;
L  = size(OldTable,1);
for i = 1:L
    
    % if trial is center out or failed, ignore
    if contains(TT.result(i),'F') || contains(TT.task(i),'CO')
        continue
    end
    
    Z = TT.RMSRatio(i);
    if Z>Threshold_up
        TT.ControlPolicy(i) = {'Velocity'};
    elseif Z<Threshold_low
        TT.ControlPolicy(i) = {'Position'};
    elseif (Z>=Threshold_low && Z<=Threshold_up)
        TT.ControlPolicy(i) = {'Hybrid'};
    end
    
end

NewTable = TT;

%% 
figure('Name','RMS Ratio')
clf
Labels = NewTable;
ii = ~isnan(Labels.RMSRatio);
ip = contains(Labels.ControlPolicy , 'Position');
iv = contains(Labels.ControlPolicy , 'Velocity');
X  = 1:length(ii);

subplot(2,2,1)
hold all
h1=plot(Labels.CursRMS_p(ii),Labels.CursRMS_v(ii),'.k');
h2=plot(Labels.CursRMS_p(ip),Labels.CursRMS_v(ip),'.r');
h3=plot(Labels.CursRMS_p(iv),Labels.CursRMS_v(iv),'.b');
xlabel('RMS_{Pos}')
ylabel('RMS_{vel}')
legend([h1,h2,h3],'Hybrid','Position Control','Velocity Control')

subplot(2,2,2)
hold all
h1=plot(X(ii),Labels.RMSRatio(ii),'.k');
h2=plot(X(ip),Labels.RMSRatio(ip),'.r');
h3=plot(X(iv),Labels.RMSRatio(iv),'.b');
plot([1,length(ii)],Threshold_low*[1,1],'g')
plot([1,length(ii)],Threshold_up*[1,1],'g')
xlabel('All Trials')
ylabel('RMS Ratio: (RMS_p / RMS_v)')
legend([h1,h2,h3],'Hybrid','Position Control','Velocity Control')




