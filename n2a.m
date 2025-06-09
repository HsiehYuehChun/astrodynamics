function a = n2a(n)
% n2a: �p�⵹�w�����B�ʲv�]revs/day�^�������y�D�b���b a�]km�^
%
% Input:
%   n - Mean motion [revs per day]
%
% Output:
%   a - Semi-major axis [km]

% �`��
mu = 398600.4418;  % �a�y�ޤO�`�� [km^3/s^2]
seconds_per_day = 86400;

% �ˬd��J
if n <= 0
    error('�����B�ʲv n ���������� (revs/day)');
end

% �p��y�D�g���]��^
T = seconds_per_day / n;

% �ھ� Kepler �ĤT�w�߭p��b���b
a = (mu * T^2 / (4 * pi^2))^(1/3);
end