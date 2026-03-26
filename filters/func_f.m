function [z] = func_f(sym_x, sym_u, dt)
    %% f(x) - UAV Process Model
    % sym_x: 16x1 symbolic state vector
    % sym_u: 6x1 symbolic input vector [Delta_theta; Delta_v]
    % dt: sampling time
    
    xn = sym_x;
    % Inputs minus biases
    dAng = sym_u(1:3) - sym_x(11:13);
    dVel = sym_u(4:6) - sym_x(14:16);
    
    %% 1. Quaternion Update
    q = xn(1:4);
    vecLen = sqrt(sum(dAng.^2));
    % Small angle approximation for RotToQuat mapping
    % We use the identity: deltaQuat = [cos(0.5*theta); sin(0.5*theta)*u]
    % For symbolic simplicity, we can use the Taylor expansion or the full form
    dq = [cos(0.5*vecLen); dAng/vecLen * sin(0.5*vecLen)];
    % If vecLen is too small, dq is [1;0;0;0]
    
    % QuatMult: q_next = [q0*dq0 - q_v'*dq_v; q0*dq_v + dq0*q_v + cross(q_v, dq_v)]
    q_next = [q(1)*dq(1) - q(2:4)'*dq(2:4); ...
              q(1)*dq(2:4) + dq(1)*q(2:4) + cross(q(2:4), dq(2:4))];
          
    %% 2. Velocity Update (NED)
    q0 = q_next(1); q1 = q_next(2); q2 = q_next(3); q3 = q_next(4);
    Tbn = [q0^2 + q1^2 - q2^2 - q3^2, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2); ...
           2*(q1*q2 + q0*q3), q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1); ...
           2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), q0^2 - q1^2 - q2^2 + q3^2];
       
    v_ned_next = xn(5:7) + Tbn * dVel + [0; 0; 9.80665] * dt;
    
    %% 3. Position Update (NED)
    p_ned_next = xn(8:10) + 0.5 * dt * (xn(5:7) + v_ned_next);
    
    %% 4. Bias Update (Constant)
    wb_next = xn(11:13);
    ab_next = xn(14:16);
    
    z = [q_next; v_ned_next; p_ned_next; wb_next; ab_next];
end

