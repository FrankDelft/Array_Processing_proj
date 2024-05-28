function freqs = espritfreq(X, d)

    % Singular value decomposition
    [U, ~, V] = svd(X, 'econ');
    V1 = V(1:end-1, 1:d);
    V2 = V(2:end, 1:d);
    D = eig(pinv(V1) * V2);
    % Frequency estimation
    freqs = angle(D)  / (2 * pi);
end
