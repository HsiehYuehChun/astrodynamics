function Rdot = F_gravity(t, R, jd0, degree)
% 二體問題基礎上加入地球重力場模型（如 EGM2008）提供的非球形重力攝動
% F_gravity: Orbital dynamics including spherical harmonic gravity perturbations
%
% Inputs:
%   t      - Time (not used, included for ODE solver compatibility)
%   R      - State vector [x; y; z; vx; vy; vz] in km and km/s (ECI frame)
%   jd0    - Julian date (for Earth rotation angle)
%   degree - Maximum degree of gravity model (e.g., 2, 10, 50, 70)
%
% Output:
%   Rdot   - Time derivative of the state vector [vx; vy; vz; ax; ay; az]

% Constants
mu = 398600.4418;  % km^3/s^2, Earth's gravitational parameter
Re = 6378137;      % m, Earth's equatorial radius (WGS-84)

% Convert ECI position to ECEF position in meters
[GMST_s, ~] = get_gst(jd0);      % get Greenwich Mean Sidereal Time in seconds
theta = mod(GMST_s / 240, 360);  % convert to degrees (360 deg = 86400 s => 1 deg = 240 s)

% Rotation matrix from ECI to ECEF (around Z axis)
rot_ECI2ECEF = [cosd(theta) sind(theta) 0;
               -sind(theta) cosd(theta) 0;
                0           0          1];

pos_ECI_km = R(1:3);
pos_ECEF_m = rot_ECI2ECEF * (pos_ECI_km * 1000);  % convert km to meters

% Compute gravity acceleration in ECEF frame using EGM2008
[gx, gy, gz] = gravitysphericalharmonic(pos_ECEF_m, 'EGM2008', degree);
acc_ECEF_kmps2 = [gx; gy; gz] / 1000;  % convert from m/s^2 to km/s^2

% Convert gravity acceleration back to ECI frame
rot_ECEF2ECI = rot_ECI2ECEF';  % transpose of rotation matrix
acc_ECI_kmps2 = rot_ECEF2ECI * acc_ECEF_kmps2;

% Compute 2-body gravity (central) acceleration
r_km = norm(pos_ECI_km);
a_grav_2body = -mu / r_km^3 * pos_ECI_km;

% Combine central gravity + spherical harmonic perturbation
a_total = a_grav_2body + acc_ECI_kmps2;

% Output time derivatives
Rdot = [R(4:6); a_total];
end