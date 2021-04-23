filename = '18-00.01.23.175-lotr.wav';
[x , fs] = audioread(filename);
trill_start_time = 0.4;

xx = x(1:floor(trill_start_time*fs));
N = length(xx);
% K = floor(sqrt(N));
% M = floor(N/K);
K = 50; %number of segments
M = floor(N/K); % segment length


nfft = 2^9;

xxx = reshape(xx(1:M*K) , M , K);
xxx = [xxx ; zeros(nfft-M ,  K) ] ;

ff = abs(fft(xxx)).^2;

psdest = mean(ff , 2);

halflen = length(psdest)/2;
figure(1)
subplot(2,1,1)
plot(xx)
subplot(2,1,2)
plot((0:(halflen-1))*fs/length(psdest) , psdest(1:halflen))
% hold on
% plot((0:(halflen-1))*fs/length(psdest) , 10*log10(psdest(1:halflen)))
title('Noise PSD Estimate (in db)')
xlabel('Hz')
% plot(ff(:,1))
% hold off


% sound(xx , fs)

impres = ifft(sqrt(psdest));
% plot(impres)

%% Window length is chosen to be smaller to reduce aliasing effect in the time domain
winlen = floor(length(impres)/2);
if(~mod(winlen, 2)) , winlen = winlen+1; end
win = hamming(winlen);
agmonfilt = [impres(end-(winlen-1)/2:end) ; impres(1:(winlen-1)/2)] .* win;
plot(real(fft(impres)))

% figure(2)
% plot([impres(1:(winlen-1)/2) ; impres(end-(winlen-1)/2:end)])
% hold on
% plot(agmonfilt)
% hold off



%% Normalize Agmon noise to obtain the same power as white noise.
poweragmon = sum(abs(fft(agmonfilt)).^2)*fs/length(agmonfilt);
powerwhite = sum(ones(size(agmonfilt)))*fs/length(agmonfilt);
gain = sqrt(powerwhite/poweragmon);
agmonfilt = agmonfilt * gain;



%% plot filter freq response for thesis or something:
figure(1000)
faxis = linspace(0 , 1 , length(agmonfilt));
agfiltpsd = abs(fft(agmonfilt)).^2;
plot(faxis , 10*log10(agfiltpsd / max(agfiltpsd)));
title('Agmon Noise Filter Power Spectrum')
xlabel('Normalized Frequency (Cycles/Sample)');
ylabel('Intensity (DB)')




wnoise = randn(N,K);
pause(0.8)
sound(0.5*wnoise(:,1) , fs)

agnoise = filter(agmonfilt , 1 , wnoise);

mean_agnoise_freqres = mean(abs(fft(agnoise)).^2, 2);

figure(3)
subplot(2,1,1)
plot(wnoise(:,1))
hold on
plot(agnoise(:,1))
hold off
title('white and agmon synth noise')
freqaxis = (0:(length(wnoise)-1))/length(wnoise)*fs;
subplot(2,1,2)
plot(freqaxis , mean(abs(fft(wnoise)).^2 , 2)/length(freqaxis))
hold on
plot(freqaxis , mean_agnoise_freqres/length(freqaxis))
hold off

legend({'white noise' , 'agmon noise'})

pause(0.6)
sound(agnoise(:,1) , fs)

fprintf('White noise power: %d \nAgmon noise power: %d\n',...
    sum(mean(abs(fft(wnoise)).^2 , 2))/length(freqaxis) , ...
    sum(mean_agnoise_freqres)/length(freqaxis));


% %% 3 secs of white noise:
% wnoise3 = randn(3*fs,1);
% %%3 secs of agmon noise:
% agnoise3 = filter(agmonfilt , 1 , wnoise3);
% soundsc(agnoise3,fs)
% audiowrite('agnoise.wav',agnoise3/2,fs)







