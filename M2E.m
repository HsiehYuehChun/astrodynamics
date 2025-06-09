function EA = M2E(M_deg, e)
% M2E: Solve Kepler's Equation (M = E - e*sin(E)) using Newton-Raphson method
% 
% Inputs:
%   M_deg - Mean anomaly in degrees
%   e     - Eccentricity (0 <= e < 1)
%
% Output:
%   EA    - Eccentric anomaly in degrees

% Convert mean anomaly to radians
M = deg2rad(M_deg);

% Normalize M to the range [-pi, pi] for better convergence
M = mod(M + pi, 2*pi) - pi;

% Initial guess E0
if M > pi
    E = M - e;
elseif M > -pi && M < 0
    E = M - e;
else
    E = M + e;
end

% Newton-Raphson iteration
tolerance = 1e-8;
delta = 1;
while abs(delta) > tolerance
    delta = (E - e * sin(E) - M) / (1 - e * cos(E));
    E = E - delta;
end

% Convert back to degrees
EA = rad2deg(E);
end