function Rdot = F_expon(t,R,jd0)
% 包含大氣阻力的二體運動模型
% F_expon: Two-body orbital dynamics with exponential atmospheric drag
% Inputs:
%   t    - Time (not used, required for ODE solvers)
%   R    - State vector [x; y; z; vx; vy; vz] in km and km/s
%   jd0  - Julian date at initial time (for Earth rotation calculation)
%
% Output:
%   Rdot - Time derivative of state vector

% Constants
mu = 398600;                % km^3/s^2, Earth gravitational parameter
Re = 6378.137;              % km, Earth equatorial radius
Cd = 2.2;                   % Drag coefficient
A_m = 4.468e-9;             % Area-to-mass ratio in m^2/kg
omega_earth = 7.2921159e-5; % rad/s, Earth rotation rate

% Position and velocity vectors
pos_km = R(1:3);            % km
vel_kmps = R(4:6);          % km/s

% Convert to meters and m/s for drag calculations
pos_m = pos_km * 1e3;
vel_mps = vel_kmps * 1e3;

% Altitude above Earth's surface in km
r_km = norm(pos_km);
alt_km = r_km - Re;

% Exponential atmosphere model (valid for ~150–1000 km)
h0 = 700;                           % reference height [km]
rho0 = 3.614e-5;                    % density at h0 [kg/m^3]
H = 88.667;                         % scale height [km]
rho = rho0 * exp(-(alt_km - h0) / H);  % atmospheric density [kg/m^3]

% Earth's rotation vector
omega_vec = [0; 0; omega_earth];   % rad/s

% Relative velocity (ECI velocity minus Earth rotation at position)
v_rel_mps = vel_mps - cross(omega_vec, pos_m);  % m/s
v_rel_mag = norm(v_rel_mps);

% Drag acceleration in m/s^2
a_drag_mps2 = -0.5 * Cd * A_m * rho * v_rel_mag * v_rel_mps;

% Convert drag acceleration back to km/s^2
a_drag_kmps2 = a_drag_mps2 / 1e3;

% Gravitational acceleration
a_grav_kmps2 = -mu / r_km^3 * pos_km;

% Total acceleration
acc_total_kmps2 = a_grav_kmps2 + a_drag_kmps2;

% Return derivative of state vector
Rdot = [vel_kmps; acc_total_kmps2];

end