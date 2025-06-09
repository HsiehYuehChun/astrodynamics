function [R, V] = randv(a, e, i, omega, w, nu)
% randv: Convert classical orbital elements to ECI state vectors
%
% Inputs:
%   a      - Semi-major axis [km]
%   e      - Eccentricity
%   i      - Inclination [deg]
%   omega  - RAAN (Right Ascension of Ascending Node) [deg]
%   w      - Argument of Periapsis [deg]
%   nu     - True anomaly [deg]
%
% Outputs:
%   R      - Position vector in ECI [km]
%   V      - Velocity vector in ECI [km/s]

% Step 1: Get position and velocity in PQW (Perifocal) frame
[r_pqw, v_pqw] = rv_pqw(a, e, nu);  % r_pqw, v_pqw ¬° 2D ¦V¶q

% Step 2: Build rotation matrix from PQW to ECI
Rz_omega = [cosd(omega), -sind(omega), 0;
            sind(omega),  cosd(omega), 0;
                 0     ,       0     , 1];

Rx_i = [1,      0       ,       0      ;
        0, cosd(i), -sind(i);
        0, sind(i),  cosd(i)];

Rz_w = [cosd(w), -sind(w), 0;
         sind(w),  cosd(w), 0;
           0   ,    0     , 1];

% Total rotation: ECI = Rz(omega) * Rx(i) * Rz(w) * r_pqw
Q_pqw2eci = Rz_omega * Rx_i * Rz_w;

% Step 3: Convert r and v to 3D and apply rotation
r_pqw_3d = [r_pqw(1); r_pqw(2); 0];
v_pqw_3d = [v_pqw(1); v_pqw(2); 0];

R = Q_pqw2eci * r_pqw_3d;
V = Q_pqw2eci * v_pqw_3d;

end