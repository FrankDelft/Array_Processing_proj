
function [X,A,S] = gen_data(M,N,delta,theta,f,SNR)
    
    
    
    % Initialize matrices
    X = zeros(M, N);
    A = zeros(M, length(theta));
    S = zeros(length(theta), N);
    num_sources=length(theta);
    Noise=zeros(M,N);
    
    % Generate the signals and array responses
    for d = 1:num_sources
        %calculate signal amplitude
        SNR_linear = 10.^(SNR(d) / 10);
        amplitude =sqrt(SNR_linear);
        for n = 1:N
            % Generate signal
            S(d, n) = amplitude*exp(1i * 2 * pi * f(d) * n);
        end
        for m = 1:M
            % Calculate array response
            A(m, d) = exp(1i * 2 * pi * delta * sin(theta(d)));
        end
    end

    % Generate noise with consistent dimensions
    real_noise = randn(size(X));
    imag_noise = randn(size(X));
    % Combine real and imaginary parts to form complex noise
    noise = real_noise + 1i * imag_noise;
    % Calculate the power of the generated noise
    noise_power = mean(abs(noise(:)).^2);
    % Normalize the noise to have power 1
    normalized_noise = noise / sqrt(noise_power);

    Noise= normalized_noise;

    %caculate the M by N matrix X
    X=A*S;
    X = X + Noise;
end
