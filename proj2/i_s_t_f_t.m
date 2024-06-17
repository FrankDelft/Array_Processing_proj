
function signal = i_s_t_f_t(S, L, Overlap)

    [M, K, ~] = size(S);
    D = round(L * (1 - Overlap));  % Step size
    N = L + (K - 1) * D;           % Total length of the reconstructed signal
    signal = zeros(M, N);
    win = hann(L, 'periodic')';
    win_sum = zeros(1, N);         % To normalize the windowing effect
    
    for m = 1:M
        for k = 1:K
            start = (k - 1) * D + 1;
            stop = start + L - 1;
            segment = real(ifft(squeeze(S(m, k, :))).');
            signal(m, start:stop) = signal(m, start:stop) + segment ;
        end
    end
    
end