function [ f0_yinbird , dips_yinbird , time_yinbird ] = yinbird(signal , fs ,  windur , threshold , stepsize , energy_threshold , minfreq)
% an implementation of yin bird, a yin based fundamental frequency
% estimator sepcialized for bir sounds.
% for more information, please refer to:
% C.  O’Reilly,  N.  M.  Marples,  D.  J.  Kelly,  and  N.  Harte,  
% “Yin-bird:  Improved  pitchtracking for bird vocalisations.,” in
% INTERSPEECH, pp. 2641–2645, 2016.
% inputs:
% signal: input signal
% fs: sample rate 
% windur: analysis window duration in seconds
% threshold: threshold for yin dips 
% stepsize: interval between analysis points in samples
% energy_threshold: energy threhsold for prominent frequencies
% minfreq: minimal value for fundamental frequency in Herz.
%
% outputs:
% f0_yinbird : estimated fundamental frequency at analysis points
% dips_yinbird: dip height for each analysis point
% time_yinbird: times of respective analysis points

winsize = floor(windur*fs);
overlap = winsize - stepsize;

low_pass = fir1(48 , 0.5,'low');
signal = filtfilt(low_pass ,1, signal);

[s , f , t_spect] = spectrogram(signal,0.01*fs,overlap,[], fs,'yaxis');
P = s .* conj(s);
[Pf_prom , maxind] = max(P);
fprom = f(maxind);
fprom(Pf_prom<mean(Pf_prom)*energy_threshold) = NaN;
fprom(fprom < minfreq) = NaN;


segsize = 3000;
segdur = segsize/fs;
M = ceil((length(signal) - winsize) / segsize);
frames = 0;
fprom_min = nan(M , 1);

for i=1:M
    first = frames(end)+1;
    last = find(t_spect<segdur*i , 1 , 'last');
    frames = first:last;
    try
        fprom_min(i) = min(fprom(frames));
    catch ME
        M = M-1;
        fprom_min = fprom_min(1:M);
    end
end

fprom_min = 100*floor(fprom_min/100);
minfreqs = unique(fprom_min(~isnan(fprom_min)));
minfreqs = wrev(minfreqs);

f0_yin = [];
time_yinbird = [];
dips_yin = [];
for i=1:length(minfreqs)
    disp([num2str(i),'/', num2str(length(minfreqs)) , ') running yin with minfreq ' , num2str(minfreqs(i))]);
    [f0 , dips , t] = yin3(signal , fs , windur , threshold , minfreqs(i) , stepsize);
    if(length(t) > length(time_yinbird)) , time_yinbird = t; end
    lendiff = length(f0_yin) - length(f0);
    if(lendiff > 0)
        f0 = [f0, nan(size(1:lendiff))];
        dips = [dips, nan(size(1:lendiff))];
    
    elseif(i>1 && lendiff < 0)
        f0_yin = [f0_yin , nan(size(f0_yin , 1) , -lendiff)];
        dips_yin = [dips_yin, nan(size(dips_yin , 1) , -lendiff)];
    end
    
    f0_yin(i,:) = f0;
    dips_yin(i,:) = dips;
end

f0_yinbird = nan(size(f0_yin(1,:)))';
dips_yinbird = nan(size(f0_yinbird));
start_seg = 1;
for i=1:M
    end_seg = find(time_yinbird < i*segdur , 1 , 'last');
    if(isnan(fprom_min(i)))
        f0_yinbird(start_seg:end_seg) = NaN;
    else
        f0_yinbird(start_seg:end_seg) = f0_yin(minfreqs==fprom_min(i),start_seg:end_seg);
        dips_yinbird(start_seg:end_seg) =  dips_yin(minfreqs==fprom_min(i),start_seg:end_seg);
    end
    start_seg = end_seg + 1;
end

f0_yinbird = f0_yinbird';
dips_yinbird = dips_yinbird';
end

