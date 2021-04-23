function [amp , ff] = desa1a(x , fs)

dx = x(2:end) - x(1:end-1);

energy = teager(x);
energyd = teager(dx);

%using backward difference
kk = 1 - energyd./(2 * energy(1:end-1));

omega = acos( kk );
amp = sqrt( energy(1:end-1) ./ ( 1 - kk.^2 ) );

ff = omega * fs/(2*pi);
