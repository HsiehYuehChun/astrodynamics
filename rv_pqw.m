function [r_pqw, v_pqw] = rv_pqw(a, e, nu_deg)
% rv_pqw: 計算 PQW 座標系下的衛星位置與速度向量
%
% Inputs:
%   a       - Semi-major axis [km]
%   e       - Eccentricity
%   nu_deg  - True anomaly [degrees]
%
% Outputs:
%   r_pqw   - Position vector in PQW frame [km, 3x1]
%   v_pqw   - Velocity vector in PQW frame [km/s, 3x1]

% 地球標準引力常數
mu = 398600.4418; % km^3/s^2

% 真近點角轉為弧度
nu = deg2rad(nu_deg);

% 半通徑
p = a * (1 - e^2);

% 位置向量
r = p / (1 + e * cos(nu));
r_pqw = [r * cos(nu); r * sin(nu); 0];

% 速度向量
v_r = sqrt(mu / p) * (-sin(nu));
v_theta = sqrt(mu / p) * (e + cos(nu));
v_pqw = [v_r; v_theta; 0];

end