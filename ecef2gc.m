function [lat_gc, lon_gc, alt_gc] = ecef2gc(Recef)
%  ECEF -> Geocentric 轉換，假設地球為球體

x = Recef(1);
y = Recef(2);
z = Recef(3);

lat_gc = atan(z / sqrt(x^2 + y^2)) * 180/pi;
lon_gc = atan2(y, x) * 180/pi;
alt_gc = sqrt(x^2 + y^2 + z^2) - 6378.137;  % km, 地球赤道半徑

end