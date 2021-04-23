function [detect , info] = longtrill_syllable_detection(xx , fs , yin , harmonics)
%returns a binary decision vector "detect" , which is the same size as xx
% xx = recording which contains 1 long trill and possibly some background
% noise and other sound events
% fs= sample rate
% yin = data structure of YIN pitch detection algorithm output (f0 , time ,
% dips)
% harmonics = harmonics at time instances as in yin.time vector

oldfs = fs;
oldxx = xx;

if(oldfs > 22.05*10^3)
    [xx , time , fs] = downsample2(oldxx , oldfs);
end



min_segment_time = 0.018;
[periodictimeyin , nonperiodictimeyin , persegs , nonpersegs] = ...
    pertime(yin.dips , yin.time , min_segment_time , yin.f0 , [2000 , 3500] );

[istart , iend] = longtrill_bulk_of_mass(persegs , yin);

evad = energy_vad(xx , fs , 0.01 , 0.5 , 1000 , yin.time(istart:iend));
tmp = zeros(size(yin.time));
tmp(istart:iend) = evad;
evad = tmp;

hold on
plot(yin.time , periodictimeyin , 'k')
plot(yin.time , evad , 'c')
plot(yin.time , periodictimeyin & evad , 'm')
hold off

[minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2] = harmonic_boundaries(yin , harmonics , periodictimeyin , [2000 , 3500]);
fbounds = [minf0 , maxf0; minf1 , maxf1 ; minf2, maxf2];

[info.time , info.fs ,info.xxds] = deal(time  ,fs , xx);
info.periodictime = segments2logical(persegs , info.time , info.fs);
info.nonperiodictime = segments2logical(nonpersegs , info.time , info.fs);
info.fbounds = fbounds;

% [periodictimeyin , nonperiodictimeyin , persegs , nonpersegs] = ...
%     pertime(yin.dips , yin.time , min_segment_time , yin.f0 , [2000 , 3500]);
% pt = segments2logical(persegs , info.time , info.fs);
% plot(gca , time , pt, 'b' ,yin.time , evad , 'c' ,time , info.periodictime, 'g');

for filternum=1:3
    xxfilt = filter_signal(xx , fs , fbounds(filternum,1) , fbounds(filternum,2));
    s = extract_data(xxfilt , fs , info.periodictime , info.nonperiodictime);
    s.xx = xxfilt;
    s.passband = fbounds(filternum,:);
    fbands(filternum) = s;
end

info.fbands = fbands;

% c = calc_env_coeffs([fbands.envpermed] , [fbands.snr]);
c=[1,0,0];

[high , low] = detection_thresholds(fbands , fs , c , info.periodictime , yin);
weighted_bands = weighted_sum_signal([fbands.xx], c);
wenergy = hamming_energy(weighted_bands, fs ,  0.01 , 1);
wenv = calc_envelope(weighted_bands);
detect = detect_syllables(wenergy,  high, low);
detect = postprocess_detection(fs , yin, wenergy , high, low , detect);

[info.wenv , info.highthresh, info.lowthresh , info.energy] = deal(wenv ,high , low , wenergy);
info.c = c;


function vad = energy_vad(xx , fs , windur , threshold , cutoff , time)
%time : optional vector. calculate only on time instants
%should contain timestamps less than length of xx in seconds

if(nargin < 6)
    time = 0:1/fs:((length(xx)-1)/fs-windur);
    totalE = sum(xx.^2)/length(xx);
else
    xxindices = floor(time(1)*fs):floor(time(end)*fs);
    totalE = sum(xx(xxindices).^2)/length(xxindices);
    time = time - windur/2; %move all time instances from center to beginning of window
    time(time>=length(xx)/fs-windur) = nan;
    time(time<0) = nan;
end

    
winlen = floor(windur*fs);

vad = zeros(size(time));

b = fir1(56 , cutoff/fs*2 , 'high');
xx = filter(b , 1 , xx);
transient = mean(grpdelay(b,1));
xx = [xx(transient+1:end) ; zeros(transient , 1)]; 


time = 1+floor(time*fs);


for i=1:length(time)
    xxind = time(i);
    if(isnan(xxind)) , continue; end
    win = xx(xxind:xxind+winlen-1);
    winE = sum(win.^2)/winlen;
    if(winE > totalE * threshold)
        vad(i) = 1;
    end
end

    


function [dsxx , dstime , dsfs] = downsample2(xx , fs)
downsample_factor = 2;
dsfs = fs/downsample_factor;
dsxx= resample(xx , 1 , downsample_factor);
dstime = (0:length(dsxx)-1)/dsfs;

function [minf , maxf] = bw_boundaries (fvect , selection , maxbw , p)
%fvect = vector of frequencies
%selection = vector of 1's and 0's (indices chosen)
%p = (optional) precentile above minimum 

if (nargin < 4) , p=0; end

perf = fvect(selection);
minf = prctile(perf , p);
maxf = prctile(perf , 100-p);

[minf , maxf] = clip_to_min_bw(minf , maxf , maxbw);


function [minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2] = harmonic_boundaries(yin , harmonics , periodictimeyin , f0boundaries)


[minf0 , maxf0] = bw_boundaries (yin.f0 , periodictimeyin, 800, 10);
[minf1 , maxf1] = bw_boundaries (harmonics(:,1) , periodictimeyin, 800, 0);
[minf2 , maxf2] = bw_boundaries (harmonics(:,2) , periodictimeyin, 800, 0);
if(isnan(minf2)) , minf2 = 3*minf0; end
if(isnan(maxf2)) , maxf2 = 3*maxf0; end
if(maxf0 > f0boundaries(2))
    maxf0=f0boundaries(2);
    if(minf0<f0boundaries(1))
        minf0 = f0boundaries(1);
    end
end


function xxfilt = filter_signal(xx , fs , flow , fhigh)
%filtering the signal to desired spectrum band
wn = [flow , fhigh]/(fs/2);
wn(end) = min(wn(end) , 0.9);
[b,a] = butter(5 , wn , 'bandpass');
% [b , ~] = fir1(72 , wn , 'bandpass'); a=1;
xxfilt = filtfilt(b , a , xx);

function data = extract_data(xx , fs , periodictime , nonperiodictime)
    
%calculating envelope
env = calc_envelope(xx);
%median of periodic envelope
menvper = env;
menvper(~periodictime) = nan;
data.menvper = medfilt1(menvper , floor(0.20*fs) , 'omitnan');

envpermed = median(env(periodictime)); %average periodic energy
envnonpermed = median(env(nonperiodictime)); %average nonperiodic energy
%SNR estimate
data.snr= 20*log10(envpermed/envnonpermed-1);
data.env = env;
data.envpermed = envpermed;
data.envnonpermed = envnonpermed;
data.energy = hamming_energy(xx , fs ,  0.01 , 1);

function env = calc_envelope(xx)

% env = envelope(xx);
env = envelope(xx , 200, 'peaks');

%changed on 3/11/17
% b = fir1(156 , 0.05);
% delay = mean(grpdelay(b));

% env = filter(b , 1 , env);

% env = medfilt1(env , 90);

% b = fir1(340 , 0.01);
% delay = delay + mean(grpdelay(b));

% env = filter(b , 1 , env);

b = fir1(156 , 0.02);
delay = mean(grpdelay(b));
env = filter(b , 1 , env);

% figure(100)
% freqz(b,1)

env = [env(delay+1:end); zeros(delay,1)];    


function [high , low] = detection_thresholds(fbands , fs , c , periodictime, yin)

energy = weighted_sum_signal([fbands.energy], c);

movmedenvper = energy;
movmedenvper(~periodictime) = nan;

[~ , badprecent] = mean_yinprobs_trill(yin);
windur = 0.02;
if(badprecent > 0.02)
    high = medfilt1(movmedenvper , floor(0.2*fs) , 'omitnan');
    low = prcfilt(energy , floor(windur*4*fs) , 40 , 'omitnan');
else
    %tested for "fast long trills" , or trills with high lobes between
    %syllables
    disp(['fast trill , ' , num2str(badprecent)])
    high = prcfilt(movmedenvper , floor(0.20*fs) , 75 , 'omitnan');
%     low = prcfilt(energy , floor(windur*4*fs) , 40 , 'omitnan');
    
%     movmedenvper = prcfilt(movmedenvper , floor(0.10*fs) , 85 , 'omitnan');
    low = medfilt1(energy , floor(windur*4*fs) , 'omitnan');
end
high(high==0) = nan;








function detect = longtrill_final_detect(fbands , fs , yin , c , periodictime)
%xx_env: cell with 3 spots , each holding envelope of f0 and the first 2
%harmonics
% envpermed: 3-array containing median of envelope for periodic parts of
% signal
%envnonpermed: 3-array containing median of envelope for non-periodic parts
%of signal
%SNR_est: SNR estimation of each of the envelopes' bandwidths
%periodictime: vector of logicals for periodic samples
%fs : sample rate
%c: optional argument, weight coefficients for each envelope. if not
%supplied they are calculated inside the function
xx_env = [fbands.env];
xx_filter = [fbands.xx];

wenv= weighted_sum_signal(xx_env , c);





[detect , low] = syllable_detection(weighted_bands , wenv, low , high ,periodictime , fs , energy);
% detect = cut_env_tails(weighted_bands , fs, energy , detect);
 detect = syllable_detection_postproc(weighted_bands , fs , yin, energy ,high, low , detect);
 
 
 function env_out = weighted_sum_signal(xxx, coeffs)
env_out = coeffs(1) * xxx(:,1) + coeffs(2) * xxx(:,2) + coeffs(3) * xxx(:,3);

function c = calc_env_coeffs(emedper , snr_est)

for i=1:length(snr_est)
    if(isreal(snr_est(i))) , snr_est(i) = real(snr_est(i)); 
    else snr_est(i) = 0;
    end
end
snr_est(snr_est<0) = 0;
snr_factor = snr_est;
snr_factor(snr_factor<2) = 0;
%changed on 6/11
% snr_factor(snr_factor~=0) = 1;
snr_factor = snr_factor/snr_factor(1); %snr of 30 or above is cosidered fully accurate
snr_factor

emedper = emedper(1)./emedper;
c = emedper .* snr_factor;

function energy = hamming_energy(xx , fs ,  windur , step)
energy = zeros(size(xx));
winsize = floor(windur*fs);
win = hamming(winsize);
xx = xx(:);
L = length(xx);
transient = ceil(winsize/2);
xx = [zeros(transient , 1) ; xx ; zeros(transient , 1)];
for i=1:L
    xxwin = xx(i:i+winsize-1);
    e = win .* xxwin;
    e = sum(e.^2);
    energy(i) = e;
end
% energy = energy/winsize;



function detect = postprocess_detection(fs , yin, energy , high, low , detect)
%look for undetected syllables in the middle of trill
detectsegs = logical2segments(detect, fs);
seglen = detectsegs(2,:) - detectsegs(1,:);
rhythm = detectsegs(1,2:end) - detectsegs(1,1:end-1);
meanr = mean(rhythm);
stdr = std(rhythm);
gaps = abs(rhythm-meanr)> stdr*2;
for i=1:length(gaps)
    if(gaps(i)>0)
        t = (detectsegs(1,i)+detectsegs(1,i+1))*0.5;
        meanl = mean(seglen(i:i+1));
        newseg = [t , t+meanl];
        newsegyin = yin.time>newseg(1) & yin.time<newseg(2);
        if(mean(yin.dips(newsegyin)) < 0.2 && std(yin.dips(newsegyin)) < 0.15)
            disp(['new syllable spotted:' , num2str(newseg(1)) , ':', num2str(newseg(2))]);
            indices = floor([detectsegs(2,i),detectsegs(1,i+1)].*fs);
            indices = indices(1):indices(2);
            detectsyl = detect_syllables(energy(indices) , high(indices)/4 , low(indices));
            detect(indices) = detectsyl;
        end
    end
    
end


%look for undetected syllables at the end of trill
tryflag = true;
while(tryflag)
    t = detectsegs(1,end)+rhythm(end); %approximate start of last missing syllable
    newseg = [t , t+seglen(end)];
    newsegyin = yin.time>newseg(1) & yin.time<newseg(2);
    if(mean(yin.dips(newsegyin)) < 0.3 && std(yin.dips(newsegyin)) < 0.3)
        disp(['new syllable spotted:' , num2str(newseg(1)) , ':', num2str(newseg(2))]);
        indices = floor([detectsegs(2,end),detectsegs(2,end)+meanr+2*stdr].*fs);
        indices(2) = min(indices(2) , length(energy));
        indices = indices(1):indices(2);
        highnan = find(isnan(high(indices)));
        if(~isempty(highnan)) %if high(indices) contains nan
            high(indices(highnan)) = high(indices(highnan(1))-1);
        end
            
        detectsyl = detect_syllables(energy(indices) , high(indices)/4 , low(indices));
        detect(indices) = detectsyl;
        syltime = indices(detectsyl>0);
        if(isempty(syltime)) , tryflag = false; 
        else
            syltime = [syltime(1) ; syltime(end)]/fs;
            detectsegs = [detectsegs , syltime];
        end
    else
        tryflag = false;
    end
end


function detect = detect_syllables(data ,  highthresh , lowthresh)
if(size(data) ~= size(highthresh))
    highthresh = highthresh';
end

abovemed = data > highthresh;

[peaks , ipeaks] = extrema(data , 1 , 0);

ipeaks = ipeaks(abovemed(ipeaks));


detect = zeros(size(data));


for i=1:length(ipeaks)
    istart=ipeaks(i);
    iend = ipeaks(i);
    while(istart>1 && data(istart)>=lowthresh(istart))
        istart=istart-1;
    end
    while(iend<length(data) && data(iend)>=lowthresh(iend))
        iend = iend+1;
    end
    
    
    if(lowthresh(istart) > lowthresh(iend))
        tightbound = lowthresh(istart);
        while(data(iend) < tightbound)
            iend = iend-1;
        end
    else
        tightbound = lowthresh(iend);
        while(istart<=length(data) && data(istart) < tightbound)
            istart = istart+1;
        end
    end
        
    if(istart~=iend)
        detect(istart:iend) = 1;
    end
end