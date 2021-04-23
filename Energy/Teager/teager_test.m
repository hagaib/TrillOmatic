T = 2;
fs = 10^3;
f = 60;

f0s = [f , 2*f , f];
times = [0 , T/2 , T];
amps = [1.5 , 0.5 , 1];
harmonics_amp_ratio = inf;
init_phase_setting = 0;

yy = generate_chirp(fs , f0s , times , amps , ...
    harmonics_amp_ratio , init_phase_setting);

% t = 0:1/fs:T-1/fs;
% yy = cos(2*pi*f*t);

spectrogram(yy , 256 , 128 , 2048 , fs , 'yaxis')
% sound(yy, fs)


tyy = teager(yy);

[amp , freq] = desa2(yy , fs);


figure(3)
p1 = subplot(3 , 1 , 1);
plot(yy)
title('sine');
p2 = subplot(3 , 1 , 2);
% plot(asin(sqrt(tyy))*fs/2/pi)
plot(amp);
title('amp est')
p3 = subplot(3 , 1 , 3);
plot(freq);
title('freq est');
linkaxes([p1, p2 , p3] , 'x')


[yy2 , fs2] = audioread('5-00.00.50.191-lotr.wav');
b = fir1(30 , [2000 3500]/(fs2/2) , 'bandpass');
yy2bp = filter(b , 1 , yy2);
delay = mean(grpdelay(b,1));
yy2bp = [yy2bp(delay+1:end) ; zeros(delay , 1)];
tyy2 = teager(yy2bp);

figure(4)
p1 = subplot(3,1,1);
plot(yy2)
title('5-00.00.50.191-lotr.wav');
p2 = subplot(3,1,2);
plot(yy2bp)
title('bandpass filtered bandwidth 2-3.5 khz');
p3 = subplot(3,1,3);
plot(tyy2)
title('teager energy and hilber envelope');

linkaxes([p1, p2 , p3] , 'x')