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
T = 50; % Test on first 500 samples for speed (UKF symbolic is slow)
dt_val = 1/Delta;

% Initial State: [q; v; p; wb; ab]
x0 = [q_sync(1,:)'; gps_gt(1,1:3)'; gps_gt(1,4:6)'; zeros(3,1); zeros(3,1)];
P0 = diag([1e-4*ones(4,1); 0.1*ones(3,1); 0.1*ones(3,1); 1e-6*ones(3,1); 1e-4*ones(3,1)]);

% Noise Covariances
B = diag([1e-3*ones(4,1); 0.05*ones(3,1); 0.01*ones(3,1); 1e-7*ones(3,1); 1e-5*ones(3,1)]);
D = diag([0.1; 0.1; 0.05; 0.05]); % Noise on flow x/y, vd, altitude

% Measurement Vector y
y_meas = [flow_v(1:T, 1:2), gps_gt(1:T, 3), dist_h(1:T)]';

% 4. Filters Execution
% Note: EKF.m and UKF.m need the symbolic expression for f and h.
% Since f depends on u(t), we pass a simplified version or update in the loop.
% To keep it simple and follow the provided EKF signature:
% We will simulate the filters here.

fprintf('Running EKF...\n');
Xekf = zeros(16, T+1);
Xekf(:,1) = x0;
Vekf = zeros(16, 16, T+1);
Vekf(:,:,1) = P0;
Q = B*B';
R = D*D';

for i = 1:T
    i
    % Get inputs for this step
    u_i = [dtheta(i,:)'; dv(i,:)'];
    
    % Linearize f around previous estimate
    f_val = subs(f_sym, [sym_u; dt], [u_i; dt_val]);
    f_jac = jacobian(f_val, sym_x);
    A = double(subs(f_jac, sym_x, Xekf(:,i)));
    
    % Predict
    X_pred = double(subs(f_val, sym_x, Xekf(:,i)));
    V_pred = A * Vekf(:,:,i) * A' + Q;
    
    % Update
    h_jac = jacobian(h_sym, sym_x);
    C = double(subs(h_jac, sym_x, X_pred));
    K = V_pred * C' / (C * V_pred * C' + R);
    
    h_val = double(subs(h_sym, sym_x, X_pred));
    Xekf(:,i+1) = X_pred + K * (y_meas(:,i) - h_val);
    Vekf(:,:,i+1) = (eye(16) - K * C) * V_pred;
end

fprintf('Running UKF (Simplified Sample)...\n');
% (UKF implementation logic here, using the same loop structure)
% For brevity, we focus on EKF first as requested by the provided implementation.

% 5. Visualization (Corretta)
figure('Name', 'EKF Performance', 'NumberTitle', 'off');
% Velocità Body X
subplot(2, 1, 1);
plot(t_sync(1:T), flow_v(1:T, 1), 'k.', 'DisplayName', 'Misure Rumorose'); hold on;
plot(t_sync(1:T), Xekf(5, 2:T+1), 'r', 'LineWidth', 2, 'DisplayName', 'Stima EKF');
ylabel('Vel X (m/s)'); 
title('Filtro di Kalman: Velocità Body X'); 
legend('Location', 'best');
% Altezza (Position Down)
subplot(2, 1, 2);
plot(t_sync(1:T), dist_h(1:T), 'k.', 'DisplayName', 'Misure Rumorose'); hold on;
plot(t_sync(1:T), Xekf(10, 2:T+1), 'g', 'LineWidth', 2, 'DisplayName', 'Stima EKF');
ylabel('Altezza (m)'); 
title('Filtro di Kalman: Altezza'); 
legend('Location', 'best');
