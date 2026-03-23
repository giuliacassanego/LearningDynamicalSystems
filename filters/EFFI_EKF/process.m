function [x]=process(xn,Delta_theta_m,Delta_v_m,Delta_theta_n,Delta_v_n,Delta_t)
% xn: State filter $\hat x_{t|t}$ at time t
% Delta_theta_m,Delta_v_m: Input
% Delta_theta_n,Delta_v_n: Parameters w.r.t. process noise
% Delta_t: Sampling time

x = zeros(length(xn),1);
%% Quaternion 
dAngTruth = Delta_theta_m' - xn(11:13) - Delta_theta_n;
deltaQuat = RotToQuat(dAngTruth);%Angle increment converted to quaternion
x(1:4)  = QuatMult(xn(1:4),deltaQuat);%Update quaternions
x(1:4)  = NormQuat(x(1:4));%Quaternion normalization

%% Rotation matrix
[Tbn] = Quat2Tbn(x(1:4));%Rotation matrix

%% Velocity
delVel    = Delta_v_m' - xn(14:16) - Delta_v_n;
delVelNav = Tbn * delVel + [0;0;0.98665]*Delta_t;
prevVel   = xn(5:7);
x(5:7)    = xn(5:7) + delVelNav;

%% Position
x(8:10) = xn(8:10) + 0.5 * Delta_t * (prevVel + x(5:7));

%% Offset
x(11:13) = xn(11:13);
x(14:16) = xn(14:16);


