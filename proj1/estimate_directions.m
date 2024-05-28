clear all;
%implementation of esprit
M = 5; % Number of antennas
N = 20; % Number of samples
theta = [-20 * (pi/180); 30 * (pi/180)]; % Directions of arrival in radians
f = [0.1; 0.3]; % Normalized frequencies
SNR = [100; 20]; % Signal-to-noise ratios in dB
delta=0.5
% Generate data and perform SVD for original data
[X, ~, ~] = gen_data(M, N, delta, theta, f, SNR);

theta=esprit(X,2);

theta_deg=rad2deg(asin(log(theta)/(2*pi*j*delta)));

