
clear all;
%Speech Enahncement for farfield/nearfield speech signals
%no difference in damping for each microphone, only phase difference

% Read the audio files
num_samples=500000;
[clean_speech, fs] = audioread('clean_speech.wav');
start_noise_seconds=2;
noise_samples=start_noise_seconds*fs;
clean_speech = clean_speech(1:num_samples);
clean_speech=[zeros(noise_samples,1);clean_speech];

[babble, ~] = audioread('babble_noise.wav');
babble = babble(1:num_samples);
babble=[zeros(noise_samples,1);babble];

[clean_speech_2, ~] = audioread('clean_speech_2.wav');
clean_speech_2 = clean_speech_2(1:num_samples);
clean_speech_2=[zeros(noise_samples,1);clean_speech_2];

[artificial_nonstat_noise, ~] = audioread('aritificial_nonstat_noise.wav');
artificial_nonstat_noise = artificial_nonstat_noise(1:num_samples);
artificial_nonstat_noise=[zeros(noise_samples,1);artificial_nonstat_noise];
num_samples=num_samples+noise_samples;

M=4; % Number of microphones


sigma = 0.1; % Variance of the white noise
%read in impulse responses
H = load('impulse_responses.mat');
h_target=getImpulse(H,5);
h_1=getImpulse(H,1);
h_2=getImpulse(H,2);


%construct multiple sound_sources
microphones = zeros(4, num_samples+399);
clean_mic = zeros(4, num_samples+399); 
interference_mic=zeros(4, num_samples+399);

for i = 1:4
    target_signal=conv(clean_speech,h_target(i,:));
    interferer1=conv(clean_speech_2,h_1(i,:));
    interferer2=conv(babble,h_2(i,:));
    white_noise = sigma * randn(size(target_signal));
    signal=target_signal+interferer1+interferer2+white_noise;

    microphones(i,:) = signal;
    clean_mic(i,:) = target_signal;
    interference_mic(i,:)=interferer1+interferer2+white_noise;
    
    filename = sprintf('microphone%d.wav', i);
    audiowrite(filename, microphones(i,:), fs);
end

% Short time fourier transform
window_time=32*10^-3;
L = 2^nextpow2(window_time*fs);
Overlapp=0.50;
[mic_stft,K] = s_t_f_t(microphones,L,Overlapp);
[inter_stft,~] = s_t_f_t(interference_mic,L,Overlapp);
alpha=0.1;
%calculate the R_x
R_x = zeros(K, L, 4, 4);
for k = 2:K
    for l = 2:L
        R_x(k, l, :, :) = alpha * squeeze(R_x(k, l-1, :, :)) + squeeze(mic_stft(:, k, l) * mic_stft(:, k, l)') / M;
    end
end


function H = getImpulse(impulse_response,index)    
    % Access the element from the struct
    impulses = reshape(struct2array(impulse_response), 5, 4, 400);
    H = squeeze(impulses(index,:,:));
end

%short time fourier transform
function [S,K] = s_t_f_t(signal,L,Overlap)
    % percentage overlap
    M=size(signal,1);
    N=length(signal);
    D=L*(1-Overlap);
    K=floor((N-L+D)/D);
    S=zeros(M,K,L);
    for m=1:M
        
        segmented_data = zeros(K, L);
        for i = 1:K
            segmented_data(i, :) = signal(1 + (i - 1) * (D) : (i-1) * (D)+L);
            S(m,i,:) = fft(segmented_data(i,:));
        end
    end

end