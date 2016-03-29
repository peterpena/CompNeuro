%% Kalman Filter

%% Belief Update

X = load('eegsample.mat');


R = size(X.EEG2);


% initial guess
S = 40;
P = 1;
RS = zeros(R);
for i=1:R
  RS(i) = S;
  [S, P] = KalmanFilter(S, P, X.EEG2(i)); 
end
y = 0;
plot(X.EEG2); hold on;
plot(RS); hold on;
%plot(1:1:R, y);