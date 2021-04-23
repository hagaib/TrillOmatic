% [xx , fs] = audioread('226-00.10.01.075-lotr.wav');
% [xx , fs] = audioread('17-00.01.26.973-lotr.wav');

[xx , fs] = audioread('5-00.00.50.191-lotr.wav');

[t , A , f0] = get_chirp_parameters(xx , fs , 0.01 , 0.001, 0.5, 2000);
yy = generate_chirp(fs , f0 , t , A , [30 ,45]);


windur = 0.001;
threshold = 0.1;
minfreq = 1500;
stepdur = 0.0015;

[ f0_yin , dips , time ] = yin3(xx , fs , windur , threshold , minfreq, floor(fs*stepdur));


 bpfilt = designfilt('bandpassfir','FilterOrder',48, ...
         'CutoffFrequency1',2000,'CutoffFrequency2',fs*0.45, ...
         'SampleRate',fs);
xx = filtfilt(bpfilt , xx);



figure(1)
plot(xx(:,1))
hold on
plot(yy)
hold off

figure(2)
spectrogram(yy,floor(0.004*fs),floor(0.003*fs),2048,fs,'yaxis');
hold on
plot(t , f0/1000)
hold off

sound(yy ,fs)

f0_plot = f0_yin;
figure(7)
spectrogram(xx,floor(0.004*fs),floor(0.003*fs),2048,fs,'yaxis');
hold on
plot(time , f0_plot/1000 , 'r');
f0_plot(dips>0.4) = nan;
plot(time , f0_plot/1000 , 'y');
f0_plot(dips>0.2) = nan;
plot(time , f0_plot/1000 , 'g');

hold off
