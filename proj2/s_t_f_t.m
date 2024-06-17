%short time fourier transform
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


