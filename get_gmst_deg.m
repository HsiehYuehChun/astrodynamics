function GMST_deg = get_gmst_deg(jd)
% get_gmst_deg: 計算指定 Julian Date 對應的 GMST（格林威治平恆星時）
%               以角度表示，範圍 [0, 360) 度
%
% Input:
%   jd         - Julian Date
%
% Output:
%   GMST_deg   - Greenwich Mean Sidereal Time [degrees]

% 日數（自 J2000.0 起）
D = jd - 2451545.0;

% GMST (in hours)，使用 IAU 近似公式
GMST_h = 18.697374558 + 24.06570982441908 * D;

% 取模到 0–24 小時之間
GMST_h_mod = mod(GMST_h, 24);

% 轉換為角度（1 小時 = 15 度）
GMST_deg = GMST_h_mod * 15;

end