
%over 1000 iterations
%compare the three algorithms esprit, espritfreq and joint
%Compare them in terms of mean values and standard deviations
%do this for SNR values varying from 0,4,8,...20dB, d=2 sources, M=3, N=20, theta=[-20,30],f=[0.1,0.12] 

% Parameters
SNR_values = 0:4:20;
d = 2;
M = 3;
N = 20;
theta = [-20, 30];
f = [0.1, 0.12];
tot_iter = 1000;
delta=0.5;
% Initialize results arrays
values_esprit = zeros(tot_iter,length(SNR_values), d);
values_espritfreq = zeros(tot_iter,length(SNR_values), d);
values_joint_freq = zeros(tot_iter,length(SNR_values), d);
values_joint_theta = zeros(tot_iter,length(SNR_values), d);



% Perform iterations
for i = 1:tot_iter
    for snr_idx = 1:length(SNR_values)
        snr = SNR_values(snr_idx);
        SNR=ones(1,d)*snr;
        % Generate data
        [X, A, S] = gen_data(M, N, delta, theta, f, SNR);
        theta_esprit = esprit(X, d,delta);
        values_esprit(i,snr_idx, :) = theta_esprit;
        freqs_espritfreq = espritfreq(X, d);
        values_espritfreq(i,snr_idx, :) = freqs_espritfreq;

        [theta_joint, f_joint] = joint(X, d, 12);
        values_joint_theta(i,snr_idx, :) = theta_joint;
        values_joint_freq(i,snr_idx, :) = f_joint;

    end
end

% Calculate means and standard deviations
mean_esprit = mean(values_esprit, 1);
std_esprit = std(values_esprit, 1);
mean_espritfreq = mean(values_espritfreq, 1);
std_espritfreq = std(values_espritfreq, 1);
mean_joint_theta = mean(values_joint_theta, 1);
std_joint_theta = std(values_joint_theta, 1);
mean_joint_freq = mean(values_joint_freq, 1);
std_joint_freq = std(values_joint_freq, 1);

% Plot means of Esprit and Joint Theta
figure;
plot(SNR_values, mean_esprit(:, :, 1), 'r-o', 'LineWidth', 2);
hold on;
plot(SNR_values, mean_joint_theta(:, :, 1), 'b-o', 'LineWidth', 2);
hold off;
xlabel('SNR (dB)');
ylabel('Mean Estimate');
title('Mean Estimates of Esprit and Joint Theta for theta 1');
legend('Esprit', 'Joint Theta');

% Plot means of Esprit and Joint Theta
figure;
plot(SNR_values, mean_esprit(:, :, 2), 'r-o', 'LineWidth', 2);
hold on;
plot(SNR_values, mean_joint_theta(:, :, 2), 'b-o', 'LineWidth', 2);
hold off;
xlabel('SNR (dB)');
ylabel('Mean Estimate');
title('Mean Estimates of Esprit and Joint Theta for theta 2');
legend('Esprit', 'Joint Theta');

% Plot means of EspritFreq and Joint Freq
figure;
plot(SNR_values, mean_espritfreq(:, :, 1), 'r-o', 'LineWidth', 2);
hold on;
plot(SNR_values, mean_joint_freq(:, :, 2), 'b-o', 'LineWidth', 2);
hold off;
xlabel('SNR (dB)');
ylabel('Mean Estimate');
title('Mean Estimates of EspritFreq and Joint Freq for f1');
legend('EspritFreq', 'Joint Freq');

figure;
plot(SNR_values, mean_espritfreq(:, :, 2), 'r-o', 'LineWidth', 2);
hold on;
plot(SNR_values, mean_joint_freq(:, :, 1), 'b-o', 'LineWidth', 2);
hold off;
xlabel('SNR (dB)');
ylabel('Mean Estimate');
title('Mean Estimates of EspritFreq and Joint Freq for f2');
legend('EspritFreq', 'Joint Freq');



% Plot standard deviations of Esprit and Joint Theta
figure;
plot(SNR_values, std_esprit(:, :, 1), 'r-o', 'LineWidth', 2);
hold on;
plot(SNR_values, std_joint_theta(:, :, 1), 'b-o', 'LineWidth', 2);
hold off;
xlabel('SNR (dB)');
ylabel('Standard Deviation');
title('Standard Deviations of Esprit and Joint Theta for theta 1');
legend('Esprit', 'Joint Theta');

% Plot standard deviations of Esprit and Joint Theta
figure;
plot(SNR_values, std_esprit(:, :, 2), 'r-o', 'LineWidth', 2);
hold on;
plot(SNR_values, std_joint_theta(:, :, 2), 'b-o', 'LineWidth', 2);
hold off;
xlabel('SNR (dB)');
ylabel('Standard Deviation');
title('Standard Deviations of Esprit and Joint Theta for theta 2');
legend('Esprit', 'Joint Theta');

% Plot standard deviations of EspritFreq and Joint Freq
figure;
plot(SNR_values, std_espritfreq(:, :, 1), 'r-o', 'LineWidth', 2);
hold on;
plot(SNR_values, std_joint_freq(:, :, 2), 'b-o', 'LineWidth', 2);
hold off;
xlabel('SNR (dB)');
ylabel('Standard Deviation');
title('Standard Deviations of EspritFreq and Joint Freq for f1');
legend('EspritFreq', 'Joint Freq');

figure;
plot(SNR_values, std_espritfreq(:, :, 2), 'r-o', 'LineWidth', 2);
hold on;
plot(SNR_values, std_joint_freq(:, :, 1), 'b-o', 'LineWidth', 2);
hold off;
xlabel('SNR (dB)');
ylabel('Standard Deviation');
title('Standard Deviations of EspritFreq and Joint Freq for f2');
legend('EspritFreq', 'Joint Freq');


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
            S(d, n) = amplitude*exp(1i * 2 * pi * f(d) * n);
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


function theta = esprit(X,d,delta)
    M=size(X,1);

    %let us group first M-1 and last M-1 data records
    X_1=X(1:M-1,:);
    Y=X(2:M,:);

    Z=[X_1;Y];
    [U_z,~,~]=svd(Z);
    U_z_hat=U_z(:,1:d);
    U_x_hat=U_z_hat(1:M-1,:);
    U_y_hat=U_z_hat(M:end,:);

    temp=pinv(U_x_hat)*U_y_hat;
    [T,Theta]=eig(temp);
    theta=sort(asind(angle(diag(Theta))/(2*pi*delta)));
end

function freqs = espritfreq(X, d)

    % Singular value decomposition
    [U, ~, V] = svd(X, 'econ');
    V1 = V(1:end-1, 1:d);
    V2 = V(2:end, 1:d);
    D = eig(pinv(V1) * V2);
    % Frequency estimation
    freqs = abs(sort(angle(D)  / (2 * pi), 'desc'));
end

function [theta,f] = joint(X,d,m)
    %lets create a m-smoothed X_m matrix
    N=size(X,2);
    M=size(X,1);
    Xm=zeros([m*M,N-m+1]);
    for i=1:m
        Xm((i-1)*M+1:i*M,: )=X(:,i:N-m+i);
    end

    %now lets take the SVD of Xm
    [Um,~,~]=svd(Xm);
    Um=Um(:,1:d);

    Um_phix=Um(1:M*(m-1),:);
    Um_phiy=Um(M+1:M*m,:);
    My=pinv(Um_phix)*Um_phiy;

    Um_theta_x=[];
    Um_theta_y=[];

    for i=1:m-1
        Um_theta_x=vertcat(Um_theta_x,Um(i*(M)+1:(i+1)*M-1,:));
        Um_theta_y=vertcat(Um_theta_y,Um(i*(M)+2:(i+1)*M,:));
    end
    Mz=pinv(Um_theta_x)*Um_theta_y;
    M=[My,Mz];
    [V,D]=joint_diag(M,1.0e-8);
    [~,v]=size(D);
    Phi=D(:,1:v/2);
    Theta=D(:,v/2+1:v);
    f=sort(angle(eig(Phi))/(2*pi),"desc");
    theta=sort( asind(angle(eig(Theta))/pi));

end
