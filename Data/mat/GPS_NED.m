function [gps_ned] = GPS_NED(lat, lon, alt)
% GPS_NED Convert latitude, longitude, altitude (in deg, deg, m)
%         to NED coordinates [m] with the first point as origin.

% 输入格式检查
if nargin < 3
    error('GPS_NED requires latitude, longitude, and altitude inputs.');
end

% 确保输入为列向量
lat = lat(:);
lon = lon(:);
alt = alt(:);

% 构造 LLA 数组
gps_lla = [lat, lon, alt];   % 单位: deg, deg, m

% 参考点（第一个采样点）
ref_lla = gps_lla(1, :);

% 经纬高 -> NED 
gps_ned = lla2ned(gps_lla, ref_lla, 'flat');

gps_ned = gps_ned';

end
