% Test parameters
M = 5; % Number of antennas
N = 20; % Number of samples
theta = [-50 * (pi/180); 50 * (pi/180)]; % Directions of arrival in radians
f = [0.1; 0.3]; % Normalized frequencies
SNR = [20; 20]; % Signal-to-noise ratios in dB

% Generate data
[X, A, S] = gen_data(M, N, 0.5, theta, f, SNR);

% Perform SVD and plot singular values
[U, Singular, V] = svd(X);

Singular