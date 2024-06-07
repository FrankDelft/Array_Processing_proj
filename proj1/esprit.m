function theta = esprit(X,d)
    M=size(X,1);

    %let us group first M-1 and last M-1 data records
    X_1=X(1:M-1,:);
    Y=X(2:M,:);

    Z=[X_1;Y];
    [U_z,~,~]=svd(Z);
    U_z_hat=U_z(:,1:d);
    U_x_hat=U_z_hat(1:M-1,:);
    U_y_hat=U_z_hat(M:end,:);

    %temp=(ctranspose(U_x_hat)*U_x_hat)*ctranspose(U_x_hat)*U_y_hat;
    temp=pinv(U_x_hat)*U_y_hat;
    [T,theta]=eig(temp);
    theta=diag(theta);
end


