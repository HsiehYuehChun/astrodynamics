function Rdot = F_2body(t, R)
%二體問題微分方程表示式
% t: time (not used in 2-body model, but required for ODE solver compatibility)

mu = 398600.4418;  % km^3/s^2

pos = R(1:3);       % [x; y; z]
vel = R(4:6);       % [vx; vy; vz]
r = norm(pos);

if r == 0
    error('Division by zero: position vector is zero.');
end

acc = -mu * pos / r^3;

Rdot = [vel; acc];
end
