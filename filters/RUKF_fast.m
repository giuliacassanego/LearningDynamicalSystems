function [Xukf, Vukf, th] = RUKF_fast(x0, y_meas, P0, Q, R, f_num, h_num, U, dt_val, c)
% RUKF_FAST Robust Unscented Kalman Filter (Numeric/Fast version)
%
% Inputs:
%   x0      : Initial state (n x 1)
%   y_meas  : Measurement matrix (m x T)
%   P0      : Initial covariance (n x n)
%   Q       : Process noise covariance (n x n)
%   R       : Measurement noise covariance (m x m)
%   f_num   : Numeric function handle for f(x, u, dt)
%   h_num   : Numeric function handle for h(x)
%   U       : Input matrix (p x T), columns are inputs at time i
%   dt_val  : Sampling time
%   c       : Robustness bound parameter
%
% Outputs:
%   Xukf    : Estimated state (n x T+1)
%   Vukf    : Estimated prediction covariance (n x n x T+1)
%   th      : Robust parameters over time (1 x T)

    n = length(x0);
    T = size(y_meas, 2);
    m = size(R, 1);
    
    alpha = 0.1; 
    kapa = 3 - n; 
    beta = 2;
    lambda = alpha^2 * (n + kapa) - n;
    Wm = [lambda/(n+lambda), ones(1, 2*n)/(2*(n+lambda))];
    Wc = Wm; Wc(1) = Wc(1) + (1 - alpha^2 + beta);
    
    x_pred = zeros(n, T+1);
    x_pred(:,1) = x0;
    
    V_pred = zeros(n, n, T+1);
    V_pred(:,:,1) = P0;
    
    Xukf = zeros(n, T+1);
    Xukf(:,1) = x0;
    
    th = zeros(1, T);
    
    for i = 1:T
        u_i = U(:, i);
        
        % -----------------------------------------
        % 1. MEASUREMENT UPDATE (from prediction)
        % -----------------------------------------
        P_sqrt = chol((n + lambda) * V_pred(:,:,i))';
        sigma_x = [x_pred(:,i), x_pred(:,i) + P_sqrt, x_pred(:,i) - P_sqrt];
        
        sigma_y_pred = zeros(m, 2*n+1);
        for j = 1:2*n+1
            sigma_y_pred(:,j) = h_num(sigma_x(:,j));
        end
        
        y_pred_mean = sum(Wm .* sigma_y_pred, 2);
        Py = R; 
        Pxy = zeros(n, m);
        for j = 1:2*n+1
            diff_x = sigma_x(:,j) - x_pred(:,i);
            diff_y = sigma_y_pred(:,j) - y_pred_mean;
            Py = Py + Wc(j) * (diff_y * diff_y');
            Pxy = Pxy + Wc(j) * (diff_x * diff_y');
        end
        
        K = Pxy / Py;
        x_hat = x_pred(:,i) + K * (y_meas(:,i) - y_pred_mean);
        x_hat(1:4) = x_hat(1:4) / norm(x_hat(1:4)); % Normalize quat
        Xukf(:, i+1) = x_hat; % Save updated state for output
        
        P_hat_temp = V_pred(:,:,i) - K * Py * K';
        [U_ph, S_ph, ~] = svd(P_hat_temp);
        S_ph(S_ph < 1e-10) = 1e-10;
        P_hat = U_ph * S_ph * U_ph';
        P_hat = (P_hat + P_hat') / 2; % Enforce symmetry
        
        % -----------------------------------------
        % 2. TIME PROPAGATION (Prediction for next step)
        % -----------------------------------------
        P_sqrt_hat = chol((n + lambda) * P_hat)';
        sigma_x_hat = [x_hat, x_hat + P_sqrt_hat, x_hat - P_sqrt_hat];
        
        sigma_x_pred = zeros(n, 2*n+1);
        for j = 1:2*n+1
            sigma_x_pred(:,j) = f_num(sigma_x_hat(:,j), u_i, dt_val);
        end
        
        x_pred(:,i+1) = sum(Wm .* sigma_x_pred, 2);
        x_pred(1:4,i+1) = x_pred(1:4,i+1) / norm(x_pred(1:4,i+1));
        
        P_pred_next = Q;
        for j = 1:2*n+1
            diff_x_pred = sigma_x_pred(:,j) - x_pred(:,i+1);
            diff_x_pred_trans = diff_x_pred'; % For clarity in UKF mult
            P_pred_next = P_pred_next + Wc(j) * (diff_x_pred * diff_x_pred_trans);
        end
        P_pred_next = (P_pred_next + P_pred_next') / 2;
        
        % -----------------------------------------
        % 3. ROBUSTNESS STEP
        % -----------------------------------------
        value = 1;
        t1 = 0; 
        e = eig(P_pred_next);
        r = max(abs(e));
        t2 = (1 - 10^-5) * (r)^-1;
        
        max_iters = 100;
        iters = 0;
        while abs(value) >= 10^-7 && iters < max_iters
           tt = 0.5 * (t1 + t2);
           val_mat = eye(n) - tt * P_pred_next;
           
           if rcond(val_mat) < 1e-12
               t2 = tt;
               iters = iters + 1;
               continue;
           end
           
           value = trace(inv(val_mat) - eye(n)) + log(det(val_mat)) - c;
           if value > 0
                t2 = tt;
           else
                t1 = tt;
           end
           iters = iters + 1;
        end
        th(i) = tt;
        
        % Mathematically: (P^-1 - th*I)^-1 = P * (I - th*P)^-1
        V_pred_temp = P_pred_next / val_mat;
        
        % SVD Projection to enforce strict positive-definiteness for chol()
        [U_v, S_v, ~] = svd(V_pred_temp);
        S_v(S_v < 1e-10) = 1e-10;  % Floor eigenvalues
        V_pred(:,:,i+1) = U_v * S_v * U_v';
        V_pred(:,:,i+1) = (V_pred(:,:,i+1) + V_pred(:,:,i+1)') / 2;
    end
    
    Vukf = V_pred; % Return covariances
end
