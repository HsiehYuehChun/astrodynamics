function[a,E,i,omega,w,nu] = elorb(R,V)
h = cross(R,V);
H = norm(h);
k = [0,0,1];
n = cross(k,h);
mu = 398600.4418; % km^3/s^2
e = ( ( (norm(V)^2)-(mu/norm(R)) )*R-(dot(R,V))*V)/mu;
E = norm(e);
epsilon = (norm(V)^2)/2 - mu/norm(R);

if E ~= 1
    a = -(mu/(2*epsilon));
    %p = a*(1-E^2);
else
    %p = (H^2)/mu;
    a = inf;
end

i = acosd(h(3)/H);

if norm(n) ~= 0
    omega = acosd(n(1)/norm(n));
    if n(2) < 0
        omega = 360 - omega;
    end
else
    omega = 0;  % 或 NaN，視應用需求
end

if norm(e) ~= 0
    w = acosd(dot(n,e)/(norm(n)*norm(e)));
    if e(3)<0
        w = 360 - w;
    end

    nu = acosd(dot(e,R)/(norm(e)*norm(R)));
    if dot(R,V)<0
        nu = 360 - nu;
    end
else
    w = 0;    % 或 NaN
    nu = acosd(dot(R,V)/(norm(R)*norm(V)));  % 用方向估計
end



