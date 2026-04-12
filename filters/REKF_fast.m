function [Xrekf, Vekf, th] = REKF_fast(x0, y_meas, P0, Q, R, f_num, f_jac_num, h_num, h_jac_num, U, dt_val, c)
% REKF_FAST Robust Extended Kalman Filter (Numeric/Fast version)
%
% Inputs:
%   x0        : Initial state (n x 1)
%   y_meas    : Measurement matrix (m x T)
%   P0        : Initial covariance (n x n)
%   Q         : Process noise covariance (n x n)
%   R         : Measurement noise covariance (m x m)
%   f_num     : Numeric function handle for f(x, u, dt)
%   f_jac_num : Numeric function handle for Jacobian of f wrt x
%   h_num     : Numeric function handle for h(x)
%   h_jac_num : Numeric function handle for Jacobian of h wrt x
%   U         : Input matrix (p x T), columns are inputs at time i
%   dt_val    : Sampling time
%   c         : Robustness bound parameter
%
% Outputs:
%   Xrekf     : Estimated state (n x T+1)
%   Vekf      : Estimated covariance (n x n x T+1)
%   th        : Robust parameters over time (1 x T)

    n = length(x0);
    T = size(y_meas, 2);
    
    Xrekf = zeros(n, T+1); 
    Xrekf(:,1) = x0;
    
    Vekf = zeros(n, n, T+1); 
    Vekf(:,:,1) = P0;
    
    th = zeros(1, T);

    for i = 1:T
        i
        u_i = U(:, i);
        
        % --- PREDICTION ---
        X_pred = f_num(Xrekf(:,i), u_i, dt_val);
        A = f_jac_num(Xrekf(:,i), u_i, dt_val);
        P_pred = A * Vekf(:,:,i) * A' + Q;
        
        % Normalize Quaternions early just in case
        X_pred(1:4) = X_pred(1:4) / norm(X_pred(1:4));
        
        % --- UPDATE ---
        C = h_jac_num(X_pred);
        S = C * P_pred * C' + R;
        
        % Kalman Gain
        K = P_pred * C' / S;
        
        % Measurement Residual
        h_val = h_num(X_pred);
        y_residual = y_meas(:,i) - h_val;
        
        % State Update
        X_upd = X_pred + K * y_residual;
        
        % Normalize quaternion
        X_upd(1:4) = X_upd(1:4) / norm(X_upd(1:4));
        Xrekf(:, i+1) = X_upd;
        
        % Covariance Update (Standard)
        V_upd_temp = (eye(n) - K * C) * P_pred;
        [U_vu, S_vu, ~] = svd(V_upd_temp);
        S_vu(S_vu < 1e-10) = 1e-10;
        V_upd = U_vu * S_vu * U_vu';
        V_upd = (V_upd + V_upd') / 2; % Enforce symmetry
        
        % --- ROBUSTNESS STEP (Least-favorable covariance) ---
        value = 1;
        t1 = 0; 
        e = eig(V_upd);
        r = max(abs(e));
        t2 = (1 - 10^-5) * (r)^-1;
        
        % Bisection method to find optimal \theta
        max_iters = 100;
        iters = 0;
        while abs(value) >= 10^-7 && iters < max_iters
           tt = 0.5 * (t1 + t2);
           val_mat = eye(n) - tt * V_upd;
           
           % Safeguard against ill-conditioned matrix during search
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
        
        % Apply robust covariance inflation without double inversions
        % Mathematically: (V^-1 - th*I)^-1 = V * (I - th*V)^-1
        Vekf_temp = V_upd / val_mat; 
        
        [U_ve, S_ve, ~] = svd(Vekf_temp);
        S_ve(S_ve < 1e-10) = 1e-10;  % Floor eigenvalues
        Vekf(:,:,i+1) = U_ve * S_ve * U_ve';
        Vekf(:,:,i+1) = (Vekf(:,:,i+1) + Vekf(:,:,i+1)') / 2; % Enforce symmetry
        
    end
end
