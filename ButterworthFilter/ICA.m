%% ICA Implementation

%% Generate random data for testing
POINTS = 1000; % number of points to plot

% define the two random variables
% -------------------------------
for i=1:POINTS
	A(i) = round(rand*99)-50;              % A
	B(i) = round(rand*99)-50;              % B
end;
figure; subplot(2, 1, 1); plot(A,B, '.');                        % plot the variables
set(gca, 'xlim', [-80 80], 'ylim', [-80 80]);  % redefines limits of the graph

% mix linearly these two variables
% --------------------------------
M1 = 0.54*A - 0.84*B;                          % mixing 1
M2 = 0.42*A + 0.27*B;                          % mixing 2
subplot(2, 1, 2); plot(M1,M2, '.');                      % plot the mixing
set(gca, 'ylim', get(gca, 'xlim'));            % redefines limits of the graph

x = [M1; M2];

%% Centering Data by Subtracting Mean
mx  = mean (x');
t = mx'*ones(1, POINTS);
w_x = x-t; 

%% Whitening Data Using EVD(Eigenvalue Decomposition)

% Get the covariance matrix
x_cov = cov(w_x');

% Get the eigenvalue and eigenvector 
[E, D] = eig(x_cov);

% EVD
x_tilde = E * D^-0.5 * E' * w_x;
%figure; plot(x_tilde(1,:), x_tilde(2,:), '.');

%% Preprocess EEG Data

% center data
R = RawSignal(:, 1:100);
R = R';
[Components, Col] = size(R);
mx = mean(R');
q = mx'*ones(1, Col);
c_x = R - q;

% whiten data
c_cov = cov(c_x');
[E, D] = eig(c_cov);
R_p = E * D^-0.5 * E' * c_x;
Z = cov(R_p');

figure;
subplot(4, 1, 1);
plot(R');
title('Raw Signal')
subplot(4, 1, 2);
plot(R_p');
title('Whitened and Centered Data');

%% FastICA Algorithm
%N = Components; %number of desired components
N = Components;
C = 2;
M = Col;
X = R_p;
one_vec = ones(M, 1);
W = zeros(N, C);
for p = 1:C
    W_p = rand(N, 1);
    W_p = W_p/norm(W_p,2);
    W_0 = rand(N, 1);
    W_0 = W_0/norm(W_0,2);
   while abs(abs(W_0'*W_p)-1) > 0.000001
       %X = 2
       W_0 = W_p;
       %W_p = 1/M*X*tanh(W_p*X)'-1/M*diff(tanh(W_p*X))*one_vec*W_p;
       W_p = 1/M*X*tanh(W_p'*X)'-1/M*(1-tanh(W_p'*X).^2)*one_vec*W_p;
       for j = 1:p-1
           if p ~= 1
            W_p = W_p - W(:, j)*W_p'*W(:, j);
           end
       end
       W_p = W_p /norm(W_p, 2);
   end
   W(:, p) = W_p;
end
S =W'*X;

subplot(4,1,3);
plot(S');
title('FastICA data');

%% Filter ICA Data Using Kalman Filter
[R, C] = size(S);
RS2 = zeros(R, C);

 for i=1:R % R = 200;Q = 0.001;
     X_s = 1;
     P = 1;
    for j=1:C
          RS2(i, j) = X_s;
          [X_s, P] = KalmanFilter(X_s, P, S(i, j)); 
    end
 end

subplot(4,1,4);
plot(RS2');
title('ICA/Kalman Filtered Data');


figure;
subplot(3, 1, 1);
plot(CLEANSignal(:, 1));
title('Clean Data')
subplot(3, 1, 2);
plot(RS2(1, :));
title('ICA/Kalman Filtered Data');
subplot(3, 1, 3);
plot(RawSignal(:, 1));
title('Raw Data');