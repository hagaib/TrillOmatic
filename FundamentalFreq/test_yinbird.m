windur = 0.001;
threshold = 0.1;
stepsize = 10;
energy_threshold = 0.01;
minfreq = 2000;

filename = 'XC13683-CommonNighthawk-Chordeilesminorhenryi';

[signal,fs] = audioread([filename '.wav']);
[ f0_yinbird , dips_yinbird , time_yinbird ] = ...
    yinbird(signal , fs ,  windur , threshold , stepsize , energy_threshold , minfreq);

figure(2)
load(['Documents\sounds\Synth_birds_database\2_synth_files\sine_info\2_thrills\' filename]);
load(['Documents\sounds\Synth_birds_database\2_synth_files\GT_pitch\2_thrills\' filename]);
winsize = floor(windur*fs);
spectrogram(signal,0.01*fs,winsize-stepsize,[],fs,'yaxis');
hold on

f0_toplot = f0_yinbird;
f0_toplot(dips_yinbird>0.3) = nan;
plot(time_yinbird,f0_toplot/1000,'r');
f0_toplot = f0_yinbird;
f0_toplot(dips_yinbird>0.2) = nan;
plot(time_yinbird,f0_toplot/1000,'y');
f0_toplot = f0_yinbird;
f0_toplot(dips_yinbird>0.1) = nan;
plot(time_yinbird,f0_toplot/1000,'g');


plot(frmTime , GT_pitch/1000, 'b');
hold off