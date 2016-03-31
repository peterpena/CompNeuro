X = load('eegsample.mat');
inputData = X.EEG2;

%% Butterworth Filter
sampleRate = 500; % Hz
cutOffFreq = 1; % Hz
filterOrder = 2; % Filter order (e.g., 2 for a second-order Butterworth filter)
[b, a] = butter(filterOrder, cutOffFreq/(sampleRate/2)); % Generate filter coefficients
filteredData = filtfilt(b, a, inputData); % Apply filter to data using zero-phase filtering


%% Kalman Filter
R = size(X.EEG2);
% initial guess
S = 5;
P = 1;
RS = zeros(R);
for i=1:R
  RS(i) = S;
  [S, P] = KalmanFilter(S, P, X.EEG2(i)); 
end

plot(X.EEG2); hold on;
plot(RS); hold on;
plot(filteredData); 
xlim([0 2000]);