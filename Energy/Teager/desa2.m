function [amp , ff] = desa2(x , fs)

dx = x(3:end) - x(1:end-2);

energy = teager(x);
energy = energy(2:end-1);
energyd = teager(dx);

%using central difference 

omega = 0.5 * acos( 1 - energyd ./ (2 * energy) );
amp = 2*energy ./ sqrt(energyd);

ff = omega * fs/(2*pi);
