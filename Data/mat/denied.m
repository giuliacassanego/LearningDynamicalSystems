function [y] =  denied(y, T, I, Detla)
% T: Start time of denial
% I: Duration of the denial
% Detal: Data frequency (Here is w.r.p to IMU)

% GPS denied data 
for k=T*Detla+1:(T+I)*Detla
    y(k,:) = y(T*Detla,:);
end

end