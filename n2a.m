function a = n2a(n)
% n2a: 計算給定平均運動率（revs/day）對應的軌道半長軸 a（km）
%
% Input:
%   n - Mean motion [revs per day]
%
% Output:
%   a - Semi-major axis [km]

% 常數
mu = 398600.4418;  % 地球引力常數 [km^3/s^2]
seconds_per_day = 86400;

% 檢查輸入
if n <= 0
    error('平均運動率 n 必須為正值 (revs/day)');
end

% 計算軌道週期（秒）
T = seconds_per_day / n;

% 根據 Kepler 第三定律計算半長軸
a = (mu * T^2 / (4 * pi^2))^(1/3);
end