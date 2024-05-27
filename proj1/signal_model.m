% Test parameters
M = 5; % Number of antennas
N = 20; % Number of samples
theta = [-50 * (pi/180); 50 * (pi/180)]; % Directions of arrival in radians
f = [0.1; 0.3]; % Normalized frequencies
SNR = [20; 20]; % Signal-to-noise ratios in dB
[X, A, S] = gen_data(M, N, 0.5, theta, f, SNR);

% Perform SVD and plot singular values
[U, Singular, V] = svd(X);
[U_s, Singular_s, V_s] = svd(S);


%double samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 40; % Number of samples
[X_sample, A_sample, S_sample] = gen_data(M, N, 0.5, theta, f, SNR);

% Perform SVD and plot singular values
[ ~, Singular_sample,~ ] = svd(X_sample);


%double antennas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 10; % Number of antennas
N = 20; % Number of samples
[X_freq, A_freq, S_freq] = gen_data(M, N, 0.5, theta, f, SNR);

% Perform SVD and plot singular values
[ ~, Singular_freq,~ ] = svd(X_freq);

%make source angles closer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 5; % Number of antennas
theta = [-5 * (pi/180); 10 * (pi/180)]; % Directions of arrival in radians
[X_angle, A_angle, S_angle] = gen_data(M, N, 0.5, theta, f, SNR);

% Perform SVD and plot singular values
[ ~, Singular_angle,~ ] = svd(X_angle);



% Define the subplot arrangement
subplot_rows = 2; % Two rows for the different cases
subplot_cols = 2; % Two columns for the different cases

% Plotting the singular values
figure;

% Plot for original data
subplot(subplot_rows, subplot_cols, 1);
stem(diag(Singular));
title('Original Data');
xlabel('Index');
ylabel('Singular Value');

% Plot for doubled samples
subplot(subplot_rows, subplot_cols, 2);
stem(diag(Singular_sample));
title('Doubled Samples');
xlabel('Index');
ylabel('Singular Value');

% Plot for doubled antennas
subplot(subplot_rows, subplot_cols, 3);
stem(diag(Singular_freq));
title('Doubled Antennas');
xlabel('Index');
ylabel('Singular Value');

% Plot for closer source angles
subplot(subplot_rows, subplot_cols, 4);
stem(diag(Singular_angle));
title('Closer Source Angles');
xlabel('Index');
ylabel('Singular Value');






