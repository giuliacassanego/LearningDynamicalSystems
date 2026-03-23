function [Q] = calcQ16(wb,ab,Delta_theta_n,Delta_v_n,q0,q1,q2,q3)
% wbn: Estimate angle increment offset error
% abn: Estimate speed increment offset error
% Delta_theta_n,Delta_v_n: Parameters w.r.t. process noise
% q0,q1,q2,q3: Estimate quaternion

%% G
G=zeros(16,6);
G(1,1)=q1/2;G(1,2)=q2/2;G(1,3)=q3/2;G(2,1)=-q0/2;G(2,2)=q3/2;G(2,3)=-q2/2;
G(3,1)=-q3/2;G(3,2)=-q0/2;G(3,3)=q1/2;G(4,1)=q2/2;G(4,2)=-q1/2;G(4,3)=-q0/2;
G(5:7,4:6)=-[q0^2+q1^2-q2^2-q3^2 2*q1*q2-2*q0*q3 2*q0*q2+2*q1*q3;2*q0*q3+2*q1*q2 q0^2-q1^2+q2^2-q3^2 2*q2*q3-2*q0*q1;2*q1*q3-2*q0*q2 2*q0*q1+2*q2*q3 q0^2-q1^2-q2^2+q3^2];

%% Q
Q=G(:,:)*diag([Delta_theta_n(1)^2 Delta_theta_n(2)^2 Delta_theta_n(3)^2 Delta_v_n(1)^2 Delta_v_n(2)^2 Delta_v_n(3)^2])*G(:,:)';

%% N_process
processNoiseVariance = [zeros(1,10), (wb').^2, (ab').^2]; % zeros(1,16)].^2

kapa = 0*[zeros(1,7), 10^-10, 10^-10, 10^-10, zeros(1,6)]; % to aviod the singular matrix

for i = 1:16
    Q(i,i) = Q(i,i) + processNoiseVariance(i) + kapa(i);
end
end