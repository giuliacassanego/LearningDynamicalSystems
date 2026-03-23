function [x,P]=EKF(x,y,V,R,Delta_theta_m,Delta_v_m,Delta_theta_n,Delta_v_n,wb,ab,Delta_t,mode)
% x:       State prediction $\hat x_{t}$ at time t
% y:       Measurements $y_t$ at time t
% range:   Height from the ground
% V:       State covraiance matrix $V_t$ at time t
% R:       Measruement noise
% Delta_theta_m,Delta_v_m: Input
% Delta_theta_n,Delta_v_n: Initial parameters w.r.t. process noise
% wb,ab initial value of speed increment offset error
% Delta_t: Sampling time

%% C_t Mode selection
if strcmp(mode, 'gps')
    C = [zeros(3,4), eye(3), zeros(3,3), zeros(3,3), zeros(3,3);
        zeros(3,4), zeros(3,3), eye(3), zeros(3,3), zeros(3,3)]; %w.r.t p_NED
elseif strcmp(mode, 'flow')
    %% Linearizing optical flow measurements
    % need to figure out how to use range data: incorporate in y? put it in x?
    % just pass it as argument?
    % range = y(5);
    % y = y(1:3);
    % 
    % dh_dq0 = [(-2*x(4)*x(5) + 2*x(1)*x(6) + 2*x(2)*x(7)) * (1/range);
    %     (2*x(1)*x(5) + 2*x(4)*x(6) - 2*x(3)*x(7)) * (-1/range)];
    % 
    % dh_dq1 = [(2*x(3)*x(5) - 2*x(2)*x(6) + 2*x(1)*x(7)) * (1/range);
    %     (2*x(2)*x(5) + 2*x(3)*x(6) + 2*x(4)*x(7)) * (-1/range)];
    % 
    % dh_dq2 = [(2*x(2)*x(5) + 2*x(3)*x(6) + 2*x(4)*x(7)) * (1/range);
    %     (-2*x(3)*x(5) + 2*x(2)*x(6) - 2*x(1)*x(7)) * (-1/range)];
    % 
    % dh_dq3 = [(-2*x(1)*x(5) - 2*x(4)*x(6) + 2*x(3)*x(7)) * (1/range);
    %     (-2*x(4)*x(5) + 2*x(1)*x(6) + 2*x(2)*x(7)) * (-1/range)];
    % 
    % dh_dvn = [(2*x(2)*x(3)-2*x(1)*x(4))/range;
    %     (x(1)^2+x(2)^2-x(3)^2-x(4)^2)/(-1/range)];
    % 
    % dh_dve = [(x(1)^2-x(2)^2+x(3)^2-x(4)^2)/range;
    %     (2*x(1)*x(4)+2*x(2)*x(3))/(-1/range)];
    % 
    % dh_dvd = [(2*x(1)*x(2)+2*x(3)*x(4))/range;
    %     (2*x(2)*x(4)-2*x(1)*x(3))/(-1/range)];
    % 
    % dh_dq = [dh_dq0, dh_dq1, dh_dq2, dh_dq3];
    % dh_dv = [dh_dvn, dh_dve, dh_dvd];

    %% C_t 
    % C = [dh_dq, dh_dv, zeros(2,9);
    %     zeros(1,4), zeros(1,3), zeros(1,2), 1, zeros(1,6)]; %w.r.t p_NED
    % C = [zeros(2,4), eye(2), zeros(2,4), zeros(2,6);
    %     zeros(1,4), zeros(1,3), zeros(1,2), 1, zeros(1,6)]; %w.r.t p_NED
end

%% L_t
L = V*C'*inv(C*V*C'+R);

%% \hat x_t|t
xn= x+L*(y-C*x);

%% A and Q
[A,Q] = linearized_process(xn,...% State filter $\hat x_{t|t}$ at time t
                            Delta_theta_m,Delta_v_m,...% Input
                            Delta_theta_n,Delta_v_n,...% Initial parameters w.r.t. process noise
                            wb,ab,...% initial value of speed increment offset error
                            Delta_t... % Sampling time
                            );

%% \hat x_t+1
[x] = process(xn,...% State filter $\hat x_{t|t}$ at time t
               Delta_theta_m,Delta_v_m,...% Input
               Delta_theta_n,Delta_v_n,...% Initial parameters w.r.t. process noise
               Delta_t... % Sampling time
               );

%% V_t+1
P=A*V*A'-A*V*C'*inv(C*V*C'+R)*C*V*A'+Q; 


end


