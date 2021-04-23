fs = 44100;
freq = 2711.54;
amplitude = 1;
dur = 2;
time = 0:1/fs:dur;
phase = rand()*2*pi;
SNR = 20;
SNR_amp = 10^(SNR/20);
k=1;

fundamental = cos(2*pi*freq*time + phase);
noise = randn(size(time));
harmony1 = 0.8*rand()*cos(2*pi*freq*time*2 + rand()*2*pi);
harmony2 = 0.8*rand()*cos(2*pi*freq*time*3 + rand()*2*pi);
harmony3 = 0.8*rand()*cos(2*pi*freq*time*4 + rand()*2*pi);
signal = amplitude * (fundamental  + harmony1 + harmony2 + harmony3 + noise/SNR_amp );
signal = signal';

% figure(4)
% plot(fundamental)
% hold on
% plot(harmony1)
% plot(harmony2)
% plot(harmony3)
% plot(signal)
% hold off

% [ f0 , dips , time ] = yin3(signal , fs , 0.01 , 0.1 , 4000, 32);

figure(1)
plot(fs/length(signal)*(0:length(signal)-1) , abs(fft(signal)));


% yy1 = harmonic_plus_noise(signal , fs , freq ,0.005 , 0.0025);
yy3 = harmonic_plus_noise_bird(signal , fs ,0.005 , 0.0025);
figure(2)
plot(signal,'g');
hold on
% plot(yy1, 'b');
plot(yy3, 'r');
hold off


yy_n = signal - yy3;
figure(3)
plot(yy_n)
soundsc(yy_n,fs);