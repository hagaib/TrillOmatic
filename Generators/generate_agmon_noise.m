function noisemat = generate_agmon_noise(signal_length , cols)
% Generate Agmon Ha'Hula noise using audio file: 18-00.01.23.175-lotr.wav.
% Noise is extracted for the beginning of the recording
% output: mxn matrix where m = signal_length and n = cols
% output noise matrix is equal in power to

filter_length = 2^8;

%% load audio file
filename = '18-00.01.23.175-lotr.wav';
[x , fs] = audioread(filename);
trill_start_time = 0.4;



%% Welch's method to estimate background noise PSD 
xx = x(1:floor(trill_start_time*fs));
N = length(xx);
K = 50; %number of segments
M = floor(N/K); % segment length
nfft = 2*filter_length; % FFT size doubled to reduce aliasing effect in the time domain

xxx = reshape(xx(1:M*K) , M , K);
xxx = [xxx ; zeros(nfft-M ,  K) ] ;
ff = abs(fft(xxx)).^2;
psdest = mean(ff , 2);



%% Create Agmon Filter
impres = ifft(sqrt(psdest));

winlen = filter_length; %floor(length(impres)/2);
if(~mod(winlen, 2)) , winlen = winlen+1; end
win = hamming(winlen);
agmonfilt = [impres(end-(winlen-1)/2:end) ; impres(1:(winlen-1)/2)] .* win;
% plot(real(fft(impres)))

% Normalize Agmon Noise Filter to obtain the same power as white noise.
poweragmon = sum(abs(fft(agmonfilt)).^2)*fs/length(agmonfilt);
powerwhite = sum(ones(size(agmonfilt)))*fs/length(agmonfilt);
gain = sqrt(powerwhite/poweragmon);
agmonfilt = agmonfilt * gain;


%% Generate output
transient = length(agmonfilt);
noisemat = randn(signal_length + transient, cols);
noisemat = filter(agmonfilt , 1 , noisemat);
noisemat = noisemat(transient+1:end, :);


end