function k_indices = VAD_matlab(signal, frameLen, overlap,fs)
    % Calculate STFT parameters
    hopLen = round(frameLen * (1 - overlap)); % Hop length in samples
    K = floor((length(signal) - frameLen + hopLen) / hopLen) + 1; % Number of STFT frames
  
    % Pre-allocate k_indices for efficiency
    k_indices = zeros(K, 1);
  
    % Perform VAD using detectSpeech (replace with your implementation)
    speech_samples = detectSpeech(signal, fs)
  
    % Iterate through STFT frames
    for k = 1:K
      % Calculate frame indices
      frame_start = (k - 1) * hopLen + 1;
      frame_end = frame_start + frameLen - 1;
  
        for i = 1:length(speech_samples)
            start= speech_samples(i,1);
            endd= speech_samples(i,2);
            if (frame_start >= start && frame_end <= endd) % Check
                k_indices(k) = 1;
            end
        end
    end
  end
  