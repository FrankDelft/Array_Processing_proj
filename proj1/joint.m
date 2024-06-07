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
    f=angle(eig(Phi))/(2*pi);
    theta=rad2deg( asin(angle(eig(Theta))/pi));

end

