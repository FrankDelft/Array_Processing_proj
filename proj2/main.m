clear all;
    % Load an audio file
    [audio, fs] = audioread('microphone1.wav');

    % Normalize the signal
    audio = audio / max(abs(audio));

    % STFT parameters
    L = 320;          % Window length
    Overlap = 0.5;    % 50% overlap
    threshold = 0.005; % Threshold for VAD

    % Compute STFT
    [S, K, win] = s_t_f_t(audio', L, Overlap);  % Transpose audio to match expected dimensions

    % Perform VAD
    [~, vadResult] = vad(S, threshold);

    % Time axis for the original signal
    time = (0:length(audio)-1) / fs;

    % Time axis for the VAD result
    D = round(L * (1 - Overlap));  % Step size
    t_vad = (0:K-1) * D / fs;

    % Plotting
    figure;
% Plotting
figure;

% Plot the original signal

plot(time, audio);
xlabel('Time [s]');
ylabel('Amplitude');
legend('Audio Signal');
hold on;
% Plot the binary VAD result

plot(t_vad, (vadResult-0.5)*0.5, 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('VAD');
title('Voice Activity Detection (VAD)');



% Short-Time Fourier Transform
function [S, K, win] = s_t_f_t(signal, L, Overlap)
    M = size(signal, 1);
    N = size(signal, 2);
    D = round(L * (1 - Overlap));  % Step size
    K = floor((N - L) / D) + 1;    % Number of segments
    S = zeros(M, K, L);
    win = hann(L, 'periodic')';    % Window function
    
    for m = 1:M
        for k = 1:K
            start = (k - 1) * D + 1;
            stop = start + L - 1;
            segment = signal(m, start:stop) .* win;
            S(m, k, :) = fft(segment);
        end
    end
end

% Voice Activity Detection
function [speechBins, vadResult] = vad(S, threshold)
    % Compute energy of each frame
    energy = sum(abs(S).^2, 3);

    % Normalize energy
    energy = energy / max(energy(:));

    % VAD decision based on threshold
    vadResult = energy > threshold;

    % Initialize speech bins matrix
    speechBins = vadResult;
end
