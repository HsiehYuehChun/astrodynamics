function [GMST_s, GAST_s] = get_gst(jd)

D = jd - 2451545.0;

% Calculate GMST
GMST_h = 18.697374558 + 24.06570982441908*D;   % unit:hr
M_mod = mod(GMST_h,24);
GMST_s = M_mod*86400/24;   % unit:sec

% Calculate GAST
eps = 23.4393 - 0.0000004*D;
Omega = 125.04 - 0.052954*D;
L = 280.47 + 0.98565*D;
delta_psi = -0.000319*sind(Omega) - 0.000024*sind(2*L);
GAST_h = GMST_h + delta_psi*cosd(eps);   % unit:hr
A_mod = mod(GAST_h,24);
GAST_s = A_mod*86400/24;   % unit:sec