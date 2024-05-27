function theta = esprit(X,d)
    M=size(X,1);
    N=size(X,2);
    %let us group first M-1 and last M-1 data records
    X_1=X(1:M-1,:);
    Y=X(2:M,:);

    Z=[X_1;Y];
    %[U_z,Sigma_z,V_z]=svd(Z);
    %U_z_hat=U_z(:,1:d);

    [U_x,Sigma_x,V_x]=svd(X_1);
    [U_y,Sigma_y,V_y]=svd(Y);
    U_x_hat=U_x(:,1:d);
    U_y_hat=U_y(:,1:d);
    
    temp=(ctranspose(U_x_hat)*U_x_hat)*ctranspose(U_x_hat)*U_y_hat;
    [T,theta]=eig(temp);
    theta=diag(theta)
end

