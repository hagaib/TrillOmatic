fs = 44100;
f0s = [3000 , 2700];
times = [0 , 0.03];
amps = [0.4 0.3];
harmonics_amp_ratio = [30 50];

pulse = generate_chirp(fs , f0s , times , amps , harmonics_amp_ratio);

% spectrogram(pulse , 'yaxis');
% sound(pulse , fs)
n= 20;
pulses = cell(1,n);

for i = 1:n
    pulses{i} = pulse;
end
 
yy = generate_pulsetrain(pulses , fs , 0.07 );

sound(yy , fs)
figure(1)
plot((0:(length(yy)-1))/fs,yy)
figure(2)
spectrogram(yy,floor(0.01*fs),floor(0.0005*fs),1024,'yaxis');