fs = 44100;
freq_start = 2000.54;
freq_end = 1000.2;
amplitude = 1;
dur = 2;
time = 0:1/fs:dur;
phase = rand()*2*pi;
SNR = 20;
SNR_amp = 10^(SNR/20);

k = (freq_end - freq_start)/ dur;
ff = freq_start + k/2 * time;

% fundamental = cos(2*pi*ff.*time+phase);
% noise = randn(size(time));
% harmony1 = (0.3+rand()/2)*cos(2*pi*ff.*time*2 + rand()*2*pi);
% harmony2 = (0.3+rand()/2)*cos(2*pi*ff.*time*3 + rand()*2*pi);
% harmony3 = (0.3+rand()/2)*cos(2*pi*ff.*time*4 + rand()*2*pi);
% signal = amplitude * (fundamental  + harmony1 + harmony2 + harmony3  + noise/SNR_amp );
% signal = signal';

figure(4)
plot(abs(fft(signal(1:1000))))
% figure(4)
% plot(fundamental)
% hold on
% plot(harmony1)
% plot(harmony2)
% plot(harmony3)
% plot(signal)
% hold off

% [ f0 , dips , time ] = yin3(signal , fs , 0.01 , 0.1 , 4000, 32);


yy = harmonic_plus_noise_bird(signal , fs ,0.005 , 0.0025 , 900 , 0.1);
yy_n = signal - yy;

figure(1)
plot(signal,'g');
hold on
plot(yy, 'b');
hold off

figure(2)
% plot(noise*amplitude/SNR_amp,'g');
% hold on
plot(yy_n,'b')
% hold off

soundsc(yy_n,fs);

