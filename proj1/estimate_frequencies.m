clear all;

% Parameters
M = 5; % Number of antennas
N = 20; % Number of samples
theta = [-20 * (pi/180); 30 * (pi/180)]; % Directions of arrival in radians
f = [0.1; 0.5]; % Normalized frequencies
SNR = [100; 200]; % Signal-to-noise ratios in dB
delta = 0.5;

% Generate data
[X, ~, ~] = gen_data(M, N, delta, theta, f, SNR);

% Perform ESPRIT
freq = abs(espritfreq(X, 2));

% Display results
disp('Estimated Frequencies:');
disp(freq);
