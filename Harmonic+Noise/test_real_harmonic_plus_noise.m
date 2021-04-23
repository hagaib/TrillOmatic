% [signal,fs] = audioread('KingFisher 1007 2016-02-15 14-30 1.wav');
% [signal,fs] = audioread('XC152586-Dark-eyedJunco-Juncohyemalis.wav');
% [signal,fs] = audioread('XC44812-WhiteThroatedKingfisher08Jan2010BagaFieldsGoaIndia.mp3');
% [signal,fs] = audioread('XC19988-WTKingfisher.mp3');
% sound(signal,fs);

% [signal,fs] = audioread('5-00.00.50.191-lotr.wav');

[signal,fs] = audioread('18-00.01.49.509-lotr.wav');

signal = signal(:,1);
figure(1)
spectrogram(signal,512,256,512 ,fs ,'yaxis');
 
n = buttord(1800/fs*2 , 1400/fs*2 , 0.5 , 30); 
[z,p,g] = butter(n , 1600/fs*2 , 'high');
sos = zp2sos(z,p,g);
freqz(sos,2^10,fs)
signal = sosfilt(sos,signal);

% bl = bassline(signal);



[yy , f0 , time , harms ,  yin] = harmonic_plus_noise_bird(signal , fs ,0.001 , 0.0005  , 1800 , 0.0005);
figure(2)

plot(signal,'g');
hold on
plot(yy, 'b');
hold off
yy_n = signal - yy;

soundsc (yy,fs)
% sound (yy_n,fs)

figure(3)
spectrogram(yy,512,256,512 ,fs ,'yaxis');
hold on
plot(time , f0*10^-3,'r')
plot(yin.time , yin.f0/1000 ,'k', yin.time, harms(:,1)*10^-3 ,'k', yin.time, harms(:,2)*10^-3,'k' , yin.time, harms(:,3)*10^-3 ,'k', yin.time, harms(:,4)*10^-3,'k')
hold off
legend({'Spectrogram' , 'HNM' , 'Yin'})
