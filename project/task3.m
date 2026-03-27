%% Task 3: EKF and UKF Testing on the Optical Flow Model
%clear; clc; close all;

project_dir = fileparts(mfilename('fullpath'));
% 1. Load data
data_path = fullfile(project_dir, '..', 'Data', 'mat', 'data_sync.mat');
load(data_path);

% 2. Setup Symbolic Variables
syms q0 q1 q2 q3 vn ve vd pn pe pd wbx wby wbz abx aby abz real
syms dthx dthy dthz dvx dvy dvz dt real
sym_x = [q0; q1; q2; q3; vn; ve; vd; pn; pe; pd; wbx; wby; wbz; abx; aby; abz];
sym_u = [dthx; dthy; dthz; dvx; dvy; dvz];

% Use the functions we implemented
f_sym = func_f(sym_x, sym_u, dt);
h_sym = func_h(sym_x);

% 3. Initialization
T = 400; % Can be higher now thanks to optimization
dt_val = 1/Delta;

% Initial State: [q; v; p; wb; ab]
x0 = [q_sync(1,:)'; gps_gt(1,1:3)'; gps_gt(1,4:6)'; zeros(3,1); zeros(3,1)];
P0 = diag([1e-4*ones(4,1); 0.1*ones(3,1); 0.1*ones(3,1); 1e-6*ones(3,1); 1e-4*ones(3,1)]);

% Noise Covariances
B_mat = diag([1e-3*ones(4,1); 0.05*ones(3,1); 0.01*ones(3,1); 1e-7*ones(3,1); 1e-5*ones(3,1)]);
D_mat = diag([0.1; 0.1; 0.05; 0.05]); % Noise on flow x/y, vd, altitude

% Measurement Vector y
y_meas = [flow_v(1:T, 1:2), gps_gt(1:T, 3), dist_h(1:T)]';

% 4. Optimization: Convert Symbolic to Numeric Functions (CRITICAL for speed)
fprintf('Optimizing symbolic functions (please wait a few seconds)...\n');
f_num = matlabFunction(f_sym, 'Vars', {sym_x, sym_u, dt});
h_num = matlabFunction(h_sym, 'Vars', {sym_x});

f_jac_sym = jacobian(f_sym, sym_x);
f_jac_num = matlabFunction(f_jac_sym, 'Vars', {sym_x, sym_u, dt});

h_jac_sym = jacobian(h_sym, sym_x);
h_jac_num = matlabFunction(h_jac_sym, 'Vars', {sym_x});

% 5. Filters Execution
fprintf('Running EKF...\n');
Xekf = zeros(16, T+1); Xekf(:,1) = x0;
Vekf = zeros(16, 16, T+1); Vekf(:,:,1) = P0;
Q = B_mat*B_mat'; R = D_mat*D_mat';

tic;
for i = 1:T
    u_i = [dtheta(i,:)'; dv(i,:)'];
    % EKF Predict
    X_pred = f_num(Xekf(:,i), u_i, dt_val);
    A = f_jac_num(Xekf(:,i), u_i, dt_val);
    V_pred = A * Vekf(:,:,i) * A' + Q;
    % EKF Update
    C = h_jac_num(X_pred);
    K = V_pred * C' / (C * V_pred * C' + R);
    h_val = h_num(X_pred);
    Xekf(:,i+1) = X_pred + K * (y_meas(:,i) - h_val);
    Vekf(:,:,i+1) = (eye(16) - K * C) * V_pred;
end
t_ekf = toc;
fprintf('EKF finished in %.4f seconds.\n', t_ekf);

fprintf('Running UKF (Numerical)...\n');
Xukf = zeros(16, T+1); Xukf(:,1) = x0;
Vukf = zeros(16, 16, T+1); Vukf(:,:,1) = P0;
n = 16; m = 4;
alpha = 0.1; kapa = 3-n; beta = 2;
lambda = alpha^2 * (n + kapa) - n;
Wm = [lambda/(n+lambda), ones(1, 2*n)/(2*(n+lambda))];
Wc = Wm; Wc(1) = Wc(1) + (1 - alpha^2 + beta);

tic;
for i = 1:T
    u_i = [dtheta(i,:)'; dv(i,:)'];
    
    % 1. Sigma Points (Prediction)
    P_sqrt = chol((n + lambda) * Vukf(:,:,i))';
    sigma_x = [Xukf(:,i), Xukf(:,i) + P_sqrt, Xukf(:,i) - P_sqrt];
    
    % 2. Propagate Sigma Points (Process)
    sigma_x_pred = zeros(n, 2*n+1);
    for j = 1:2*n+1
        sigma_x_pred(:,j) = f_num(sigma_x(:,j), u_i, dt_val);
    end
    
    % 3. Predicted Mean and Covariance
    x_pred = sum(Wm .* sigma_x_pred, 2);
    P_pred = Q;
    for j = 1:2*n+1
        diff = sigma_x_pred(:,j) - x_pred;
        P_pred = P_pred + Wc(j) * (diff * diffMulti(diff)); % Custom mult for safety
    end
    P_pred = (P_pred + P_pred')/2; % Enforce symmetry

    % 4. Sigma Points (Update)
    P_sqrt_pred = chol((n + lambda) * P_pred)';
    sigma_x_upd = [x_pred, x_pred + P_sqrt_pred, x_pred - P_sqrt_pred];
    
    % 5. Propagate Sigma Points (Measurement)
    sigma_y_pred = zeros(m, 2*n+1);
    for j = 1:2*n+1
        sigma_y_pred(:,j) = h_num(sigma_x_upd(:,j));
    end
    
    % 6. Predicted Measurement and Cross-Covariance
    y_pred = sum(Wm .* sigma_y_pred, 2);
    Py = R; Pxy = zeros(n, m);
    for j = 1:2*n+1
        diff_x = sigma_x_upd(:,j) - x_pred;
        diff_y = sigma_y_pred(:,j) - y_pred;
        Py = Py + Wc(j) * (diff_y * diff_y');
        Pxy = Pxy + Wc(j) * (diff_x * diff_y');
    end
    
    % 7. Final Update
    K = Pxy / Py;
    Xukf(:,i+1) = x_pred + K * (y_meas(:,i) - y_pred);
    Vukf(:,:,i+1) = P_pred - K * Py * K';
    Vukf(:,:,i+1) = (Vukf(:,:,i+1) + Vukf(:,:,i+1)')/2; % Symmetry

    if mod(i, 100) == 0, fprintf('Step %d/%d\n', i, T); end
end
t_ukf = toc;
fprintf('UKF finished in %.4f seconds.\n', t_ukf);

function out = diffMulti(d)
    out = d';
end

% 6. Visualization
figure('Name', 'EKF vs UKF Performance', 'NumberTitle', 'off', 'Color', 'w');

% Velocità Body X
subplot(2, 1, 1);
plot(t_sync(1:T), flow_v(1:T, 1), 'k.', 'DisplayName', 'Misure Rumorose'); hold on;
plot(t_sync(1:T), Xekf(5, 2:T+1), 'r', 'LineWidth', 2, 'DisplayName', 'Stima EKF');
plot(t_sync(1:T), Xukf(5, 2:T+1), 'b--', 'LineWidth', 2, 'DisplayName', 'Stima UKF');
ylabel('Vel X (m/s)'); title('Velocity Body X: EKF vs UKF');
grid on; legend('Location', 'best');

% Altezza (Position Down)
subplot(2, 1, 2);
plot(t_sync(1:T), dist_h(1:T), 'k.', 'DisplayName', 'Misure Rumorose'); hold on;
plot(t_sync(1:T), Xekf(10, 2:T+1), 'g', 'LineWidth', 2, 'DisplayName', 'Stima EKF');
plot(t_sync(1:T), Xukf(10, 2:T+1), 'm--', 'LineWidth', 2, 'DisplayName', 'Stima UKF');
ylabel('Altezza (m)'); title('Altezza: EKF vs UKF');
grid on; legend('Location', 'best');
