function [amp , ff] = desa1(x , fs)

dx = x(2:end) - x(1:end-1);

energy = teager(x);
energy = energy(2:end-1);
energyd = teager(dx);

%using backward and forward difference average
kk = 1 - ( energyd(1:end-1) + energyd(2:end) )./(4 * energy);

omega = acos( kk );
amp = sqrt( energy ./ ( 1 - kk.^2 ) );

ff = omega * fs/(2*pi);
