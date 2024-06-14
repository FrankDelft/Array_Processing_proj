

clear all;
% Parameters
M = 5; % Number of antennas
N = 20; % Number of samples
theta = [-20 ; 30 ]; % Directions of arrival in radians
f = [0.8; 0.3]; % Normalized frequencies
SNR = [100; 200]; % Signal-to-noise ratios in dB
delta = 0.5;

% Generate data
[X, ~, ~] = gen_data(M, N, delta, theta, f, SNR);

% Perform ESPRIT
freq = espritfreq(X, 2);

% Display results
disp('Estimated Frequencies:');
disp(freq);

function freqs = espritfreq(X, d)

    % Singular value decomposition
    [~, ~, V] = svd(X);
    V1 = V(1:end-1, 1:d);
    V2 = V(2:end, 1:d);
    D = eig(pinv(V1) * V2);
    
    freqs=real(log(D)/(2*pi*j));
    for i=1:d
        if freqs(i)>0
            freqs(i)=freqs(i)-1;
        end
    end
    freqs=sort(abs(freqs));   
  
end

function [X,A,S] = gen_data(M,N,delta,theta,f,SNR)
    % Initialize matrices
    num_sources=length(theta);
    X = zeros(M, N);
    A = zeros(M, num_sources);
    S = zeros(length(theta), N);
    Noise=zeros(M,N);
    
    % Generate the signals and array responses
    for d = 1:num_sources
        %calculate signal amplitude
        SNR_linear = 10.^(SNR(d) / 10);
        amplitude =sqrt(SNR_linear);
        for n = 1:N
            % Generate signal
            S(d, n) = exp(1i * 2 * pi * f(d) * n);
        end
        for m = 1:M
            % Calculate array response
            A(m, d) = exp((m-1)*1i * 2 * pi * delta * sind(theta(d)));
        end
    end

    % Generate noise with consistent dimensions
    real_noise = randn(size(X));
    imag_noise = randn(size(X));
    % Combine real and imaginary parts to form complex noise
    Noise = real_noise + 1i * imag_noise;
    signal_power = mean(abs(S(:).^2));

    % Scale noise to achieve desired SNR for each source
    for d = 1:num_sources
        SNR_linear = 10^(SNR(d) / 10);
        noise_power = signal_power / SNR_linear;
        Noise = Noise * sqrt(noise_power / mean(abs(Noise(:)).^2));
    end
    %caculate the M by N matrix X
    X=A*S;
    X = X + Noise;
end