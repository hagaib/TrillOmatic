% this script will help debug longrtill parameters extraction
filenames = get_list_csv('filenames.csv');
[xx , fs] = audioread(filenames{2});
soundsc(xx,fs)

windur = 0.01;
stepdur = 0.001;
min_freq = 1800;
energy_threshold = 0.1;


% [yy , f0, time, harms] = harmonic_plus_noise_bird(xx , fs ,windur , stepdur , min_freq , energy_threshold);


soundsc(yy,fs)
params = extract_trill_parameters(xx , fs , f0 , time);

figure(1)
spectrogram(xx,floor(windur*fs),floor(stepdur*fs),2048,fs,'yaxis')
hold on
plot(time,f0/1000,'y')
plot(time,params.is_voiced,'r')
hold off

figure(2)
plot((0:length(xx)-1)/fs,xx)