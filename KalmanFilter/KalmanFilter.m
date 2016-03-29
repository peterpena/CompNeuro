function [ X, P ] = KalmanFilter( X_p, P_p, Z_p)
%KALMANFILTER Summary of this function goes here
%   Detailed explanation goes here

%% Initialize Matrices
R = 40;
Q = 0.000001;

%% Belief Update
X_tilde = X_p;
P_tilde = P_p + Q;

%% Measurement Update
K = P_tilde / (P_tilde + R);
X = X_tilde + K*(Z_p - X_tilde);
P = (1-K)*P_tilde;
end