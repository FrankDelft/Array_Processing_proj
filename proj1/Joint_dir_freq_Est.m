clear all;

% Parameters
M = 5; % Number of antennas
N = 20; % Number of samples
theta = [-20 * (pi/180); 30 * (pi/180)]; % Directions of arrival in radians
f = [0.1; 0.3]; % Normalized frequencies
SNR = [100; 100]; % Signal-to-noise ratios in dB
delta = 0.5;
d=2;

% Generate data
[X, ~, ~] = gen_data(M, N, delta, theta, f, SNR);



[theta_j,f_j]=joint(X,d,M)