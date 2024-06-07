clear all;
% Test parameters
M = 5; % Number of antennas
N = 20; % Number of samples
theta = [-20 * (pi/180); 30 * (pi/180)]; % Directions of arrival in radians
f = [0.1; 0.3]; % Normalized frequencies
SNR = [20; 20]; % Signal-to-noise ratios in dB

% Generate data and perform SVD for original data
[X, ~, ~] = gen_data(M, N, 0.5, theta, f, SNR);
[~, Singular] = svd(X);

% Double samples
N = 40;
[X_sample, ~, ~] = gen_data(M, N, 0.5, theta, f, SNR);
[~, Singular_sample] = svd(X_sample);

% Double antennas
M = 10; N = 20;
[X_ant, ~, ~] = gen_data(M, N, 0.5, theta, f, SNR);
[~, Singular_ant] = svd(X_ant);

% Make source angles closer
M = 5; theta = [-5 * (pi/180); 10 * (pi/180)];
[X_angle, ~, ~] = gen_data(M, N, 0.5, theta, f, SNR);
[~, Singular_angle] = svd(X_angle);

% Make frequencies closer
theta = [-20 * (pi/180); 30 * (pi/180)]; f = [0.29; 0.3];
[X_freq, ~, ~] = gen_data(M, N, 0.5, theta, f, SNR);
[~, Singular_freq] = svd(X_freq);

% Store singular values and titles in arrays
singular_values = {diag(Singular), diag(Singular_sample), diag(Singular_ant), diag(Singular_angle), diag(Singular_freq)};
titles = ["Original Data", "Doubled Samples", "Doubled Antennas", "Closer Source Angles", "Closer Frequencies"];

% Create a figure and plot all the singular values in subplots
figure;
num_plots = length(singular_values);
for i = 1:num_plots
    subplot(ceil(num_plots/2), 2, i); % Arrange subplots in a grid
    stem(singular_values{i});
    title(titles(i));
    xlabel('Index');
    ylabel('Singular Value');
end
