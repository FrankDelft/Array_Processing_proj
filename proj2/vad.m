% Voice Activity Detection
function vadResult = vad(S, threshold)
    % Compute energy of each frame
    energy = sum(abs(S).^2, 3);

    % Normalize energy
    energy = energy / max(energy(:));

    % VAD decision based on threshold
    vadResult = energy > threshold;

end
