function [x_pred, V_pred, th, x_hat, L] = RUKF(x0,y,V0,B,D,c,len)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robust UKF 
%
% Input:
%       f/h: Model
%       x_0: Initial condition x(1|0)
%       P0 : Initial condition variance P(1|0)
%       B  : Matrix B of the noise inputs
%       D  : Matrix D of the noise inputs
%       y  : Measured data
%       len: Length of the dataset, and so the length of the iterations

% Outputs:
%       x_pred: state predcited   x(t+1|t)
%       P_pred: state performance P(t+1|t)

%% SAVE
Q = B*B';
R = D*D';
n = size(V0,1);
m = size(R,1);
th= zeros(len,1);       %save \theta
x_pred = zeros(n,len+1);  %save \hat x_{t}
x_hat  = zeros(n,len);  %save \hat x_{t|t}
y_pred = zeros(m,len);  %save \hat y_{t}
V_pred = zeros(n,n,len+1); %save V_t
Wc     = zeros(1,2*n+1); %save weights Wc
Wm     = zeros(1,2*n+1); %save weights Wm
L      = zeros(n,m,len); %save Kalman gain

%% Scaling factor (see Wan and Van der Merwe, 2001)
alpha = 0.5; % small postive value (10^-4 ~ 1)
kapa  = 3-n; 
lambda= alpha^2*(kapa+n)-n;
beta  = 2;   % for Gaussian 2 is optimal

%% Initial
% First prediction (i.e. initial conditions) x(1|0)
x_pred(:,1) = x0; 

% First prediction variance P0 = P(1|0)
V_pred(:,:,1) = V0;

%% Iterative part
for t = 1:len
    %% SAVE
    X      = zeros(n,2*n+1); %save Sigma points for X^(i)_t
    Y      = zeros(m,2*n+1); %save Sigma points for Y^(i)_t=h(X^(i)_t)
    X_hat  = zeros(n,2*n+1); %save Sigma points for X^(i)_t|t
    X_pred = zeros(n,2*n+1); %save Sigma points for X_pred^(i)_t+1=f(X^(i)_t|t)
    %% Generate sigma points of x_pred_t
    for i = 1:n
        sqrtP   = chol(V_pred(:,:,t)); 
        X(:,i)  = x_pred(:,t) + sqrt(n+lambda)*sqrtP(i,:)';
        X(:,i+n)= x_pred(:,t) - sqrt(n+lambda)*sqrtP(i,:)';
    end
    X(:,2*n+1)  = x_pred(:,t); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Weights Wc/Wm
    Wc(1:2*n) = 1/(2*(n+lambda));
    Wc(2*n+1) = lambda/(n+lambda) + 1 - alpha^2 + beta;
    Wm(1:2*n) = 1/(2*(n+lambda));
    Wm(2*n+1) = lambda/(n+lambda);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Compute Y^(i)_t=h(X^(i)_t) and its mean
    for i=1:2*n+1
        Y(:,i)      = func_h(X(:,i));
        y_pred(:,t) = y_pred(:,t)+Wm(i)*Y(:,i);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Find the Kalman estimation gain
    P_y  = zeros(m,m); %sum^{2n}_0 {W^i*(Y^i-mean_Y)(Y^i-mean_Y)^T+R};
    P_xy = zeros(n,m); %sum^{2n}_0 {W^i*(X^i-mean_x)(Y^i-mean_Y)^T};
    for i=1:2*n+1
        P_y =P_y +Wc(i)*(Y(:,i)-y_pred(:,t))*(Y(:,i)-y_pred(:,t))';
        P_xy=P_xy+Wc(i)*(X(:,i)-x_pred(:,t))*(Y(:,i)-y_pred(:,t))';
    end
    P_y = P_y + R;
    L(:,:,t)   = P_xy*inv(P_y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute x(t|t)
    x_hat(:,t) = x_pred(:,t) + L(:,:,t)*(y(t)-y_pred(:,t)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute P(t|t)
    P_hat = V_pred(:,:,t) - L(:,:,t)*P_y*L(:,:,t)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate Sigma points of x_hat_t
    for i = 1:n
        sqrtP_hat   = chol(P_hat); 
        X_hat(:,i)  = x_hat(:,t) + sqrt(n+lambda)*sqrtP_hat(i,:)';
        X_hat(:,i+n)= x_hat(:,t) - sqrt(n+lambda)*sqrtP_hat(i,:)';
    end
    X_hat(:,2*n+1)  = x_hat(:,t); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute X^(i)_t+1=f(X_hat^(i)_t) and its mean
    for i=1:2*n+1
        X_pred(:,i)   = func_f(X_hat(:,i));
        x_pred(:,t+1) = x_pred(:,t+1)+Wm(i)*X_pred(:,i);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute the prediction variance P(t+1|t)
    P_pred = zeros(n,n);
    for i=1:2*n+1
        P_pred=P_pred+Wc(i)...
            *(X_pred(:,i)-x_pred(:,t+1))*(X_pred(:,i)-x_pred(:,t+1))';
    end
    P_pred = P_pred + Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% least favroable covariance matrix
    value=1;
    t1=0; 
    e = eig(P_pred);
    r = max(abs(e));
    t2=(1-10^-5)*(r)^-1;
    while abs(value)>=10^-9
       tt=0.5*(t1+t2);
       value=trace(inv(eye(n)-tt*P_pred)-eye(n)) + log(det(eye(n)-tt*P_pred))-c;
       if value>0
            t2=tt;
       else
            t1=tt;
       end
    end
    th(t) = tt;
    V_pred(:,:,t+1)=inv(inv(P_pred)-th(t)*eye(n));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



