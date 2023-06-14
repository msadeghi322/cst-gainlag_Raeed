% Calculating firing rate from spikes


% Example binary spike data
D = load('Earl_20190716_COCST_TD');
D = D.trial_data;

%%
trial = 1;
spikeData = D(trial).M1_spikes(:,1);

% Define time interval and bin size
timeInterval = 6; % in seconds
binSize = 0.05; % in seconds

% Calculate number of bins
numBins = floor(timeInterval / binSize);

% Reshape spike data into bins
spikeDataBins = reshape(spikeData(1:numBins*floor(length(spikeData)/numBins)),numBins,[]);

% Calculate firing rate for each bin
firingRate = sum(spikeDataBins, 2) / binSize;

% Plot firing rate
timeAxis = linspace(0,timeInterval,numBins);
plot(timeAxis,firingRate)
xlabel('Time (s)')
ylabel('Firing rate (spikes/s)')
title('Firing rate plot')