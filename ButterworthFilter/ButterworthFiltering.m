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

figure;
subplot(3, 2, 1);
plot(X.EEG2); hold on;
plot(RS); hold on;
plot(filteredData); 
xlim([0 2000]);
title('Butterworth/Kalman Filter');


%% Plot Data from EEGLab
subplot(3, 2, 2);
plot(RawSignal);
title('Raw Data');

subplot(3, 2, 3);
plot(CLEANSignal);
title('Clean Data');
%% Filter Raw Data Using Kalman Filter
[R, C] = size(RawSignal);
RS2 = zeros(R, C);

 for i=1:C %Works best with Q = 0.000001; R = 1;
     S = 5;
     P = 1;
    for j=1:R
          RS2(j, i) = S;
          [S, P] = KalmanFilter(S, P, RawSignal(j, i)); 
    end
 end

subplot(3, 2, 4);
plot(RS2);
title('Kalman Filtered Data')
%plot(RawSignal(:,1));

%% Filter Raw Data Using Particle Filter
C = 4;
RS3 = zeros(R, C);
 M = 25;
 for i=1:C 
    % initialize particles
    V = 100;
    X_p = ones(M, 1);
    for m = 1:M
        X_p(m) = sqrt(V) * randn;
    end
    for j=1:R
          RS3(j, i) = mean(X_p);
          X_p = ParticleFilter(X_p, RawSignal(j, i)); 
    end
 end
 
subplot(3, 2, 5);
plot(RS3(:, 1));
title('Particle Filtered Data')