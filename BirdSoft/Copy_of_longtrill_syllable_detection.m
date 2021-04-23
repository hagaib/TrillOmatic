function [detect , info] = longtrill_syllable_detection(xx , fs , yin , harmonics)
%returns a binary decision vector "detect" , which is the same size as xx
% xx = recording which contains 1 long trill and possibly some background
% noise and other sound events
% fs= sample rate
% yin = data structure of YIN pitch detection algorithm output (f0 , time ,
% dips)
% harmonics = harmonics at time instances as in yin.time vector


if(fs == 44.1*10^3)
    [xxds , dstime , new_fs] = downsample2(xx , fs);
end

min_segment_time = 0.018;
[periodictimeyin , nonperiodictimeyin , persegs , nonpersegs] = ...
    pertime(yin.dips , yin.time , min_segment_time , yin.f0 , [2000 , 3500]);

[minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2] = harmonic_boundaries(yin , harmonics , periodictimeyin , [2000 , 3500]);

[info.time , info.fs ,info.xxds] = deal(dstime  ,new_fs , xxds);
info.periodictime = segments2logical(persegs , info.time , info.fs);
info.nonperiodictime = segments2logical(nonpersegs , info.time , info.fs);
fbounds = [minf0 , maxf0; minf1 , maxf1 ; minf2, maxf2];
info.fbounds = fbounds;

for filternum=1:3
    [xx_f , env ,menvper, envpermed , envnonpermed , SNR_est , low , high] = ...
    process_filtered_signal(xxds , new_fs , fbounds(filternum,1) , fbounds(filternum,2) , info.periodictime , info.nonperiodictime);
    [s.xx , s.env , s.envpermed , s.envnonpermed , s.snr_est] = deal(xx_f , env , envpermed , envnonpermed ,SNR_est);
    s.passband = fbounds(filternum,:);
    fbands(filternum) = s;
end

info.fbands = fbands;


c = calc_env_coeffs([fbands.envpermed] , [fbands.snr_est]);
[detect , movmed , wenv , low , energy] = longtrill_final_detect(fbands, info.periodictime , new_fs , yin , c);

[info.movmed , info.wenv , info.low , info.energy] = deal(movmed , wenv , low , energy);
info.c = c;



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


function [xx_f , env , menvper , envpermed , envnonpermed , SNR_est  , low , high] = process_filtered_signal(xx , fs , lower , higher , periodictime , nonperiodictime)

%filtering the signal to desired spectrum band
wn = [lower , higher]/(fs/2);
wn(end) = min(wn(end) , 0.9);
[b,a] = butter(5 , wn , 'bandpass');
% [b , ~] = fir1(72 , wn , 'bandpass'); a=1;
xx_f = filtfilt(b , a , xx);

%calculating envelope
env = calc_envelope(xx_f);

%median of periodic envelope
menvper = env;
menvper(~periodictime) = nan;
menvper = medfilt1(menvper , floor(0.20*fs) , 'omitnan');

%SNR estimate
envpermed = median(env(periodictime)); %average periodic energy
envnonpermed = median(env(nonperiodictime)); %average nonperiodic energy
SNR_est = 20*log10(envpermed/envnonpermed-1);

%high and low bound calculations
high = envpermed;
low = envnonpermed + (envpermed - envnonpermed) * 0.2;

%median of periodic envelope
movmedenvper = env;
movmedenvper(~periodictime) = nan;
movmedenvper = medfilt1(movmedenvper , floor(0.20*fs) , 'omitnan');
movmedenvper(movmedenvper==0) = nan;

high = movmedenvper;
low = low *ones(size(env));


% energy = hamming_energy(xx_f, fs ,  0.02 , 1);
% 
% detect = syllable_detection(xx ,env , low , high , periodictime ,fs , energy);


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


function [detect , movmedenvper , wenv , low , energy] = longtrill_final_detect(fbands , periodictime , fs , yin , c)
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
envpermed = [fbands.envpermed];
envnonpermed = [fbands.envnonpermed];
SNR_est = [fbands.snr_est];
xx_filter = [fbands.xx];

if(nargin<5)
    c = calc_env_coeffs(envpermed , SNR_est);
    disp('calculating coefficients');
end

wenv= weighted_env(xx_env , c);

%high and low bound calculations
high = sum(envpermed .* c);
low = sum(envnonpermed .* c);
% p = 0.15 + 0.05*(sum(c>0));
p=0.2;
if(sum(c>0)==3)
    p=0.25;
end
low = (1-p)*low + p*high;


weighted_bands = weighted_env(xx_filter , c);
energy = hamming_energy(weighted_bands , fs ,  0.01 , 1);
movmedenvper = energy;

%median of periodic envelope
% movmedenvper = wenv;
movmedenvper(~periodictime) = nan;
[~ , badprecent] = mean_yinprobs_trill(yin);
windur = 0.02;
if(badprecent > 0.02)
    movmedenvper = medfilt1(movmedenvper , floor(0.2*fs) , 'omitnan');
    low = prcfilt(energy , floor(windur*4*fs) , 40 , 'omitnan');
else
    %tested for "fast long trills" , or trills with high lobes between
    %syllables
    disp(['fast trill , ' , num2str(badprecent)])
    movmedenvper = prcfilt(movmedenvper , floor(0.20*fs) , 75 , 'omitnan');
    low = prcfilt(energy , floor(windur*4*fs) , 40 , 'omitnan');
    
%     movmedenvper = prcfilt(movmedenvper , floor(0.10*fs) , 85 , 'omitnan');
%     low = medfilt1(energy , floor(windur*4*fs) , 'omitnan');
end
movmedenvper(movmedenvper==0) = nan;

high = movmedenvper;

% low = wenv;
% low(periodictime) = nan;
% low = medfilt1(low, floor(0.15*fs) , 'omitnan');
% low = low*ones(size(wenv));


% segs = logical2segments(periodictime , fs);
% segs_between = nan(size(segs) + [0,1]);
% segs_between(1 ,2:end) = segs(2 , 1:end) + 1/fs;
% segs_between(2 ,1:end-1) = segs(1 , 1:end) - 1/fs;
% segs_between(1,1) = 0;
% segs_between(2,end) = (length(wenv)-1)/fs;

% low = nan(size(periodictime));
% for i=1:length(segs_between)
%     between = segments2logical( segs_between(:,i),[0,(length(wenv)-1)/fs],fs);
%     low(between) = prctile(wenv(between) , 30);
% end
% 
% for i=1:length(segs)
%     pstart = floor(segs(1,i)*fs);
%     pend = ceil(segs(2,i)*fs)+1;
%     L = pend-pstart+2;
%     t = 1/L:1/L:(L-1)/L;
%     vstart = low(pstart-1);
%     vend = low(pend+1);
%     low(pstart:pend) = vstart * (1-t) + vend*t;
% end



[detect , low] = syllable_detection(weighted_bands , wenv, low , high ,periodictime , fs , energy);
% detect = cut_env_tails(weighted_bands , fs, energy , detect);
 detect = syllable_detection_postproc(weighted_bands , fs , yin, energy ,high, low , detect);
 
 
 function env_out = weighted_env(envs , coeffs)
env_out = coeffs(1) * envs(:,1) + coeffs(2) * envs(:,2) + coeffs(3) * envs(:,3);

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
c=[1,0,0];

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



function detect = syllable_detection_postproc(xx , fs , yin, energy , high, low , detect)
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
            detectsyl = syllable_detection(xx(indices) , xx(indices) , low(indices) , high(indices)/4 , [] , [] , energy(indices));
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
        indices(2) = min(indices(2) , length(xx));
        indices = indices(1):indices(2);
        highnan = find(isnan(high(indices)));
        if(~isempty(highnan)) %if high(indices) contains nan
            high(indices(highnan)) = high(indices(highnan(1))-1);
        end
            
        detectsyl = syllable_detection(xx(indices) , xx(indices) , low(indices) , high(indices)/4 , [] , [] , energy(indices));
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


function [detect , lowbound] = syllable_detection(xx , env , lowbound, high , periodictime , fs , energy)
if(size(env) ~= size(high))
    high = high';
end

abovemed = energy > high;


%peak detection
% [peaks , ipeaks] = extrema(env , 1 , 0);
[peaks , ipeaks] = extrema(energy , 1 , 0);

%give periodic peaks a detection bonus
% adjusted_high = high - (high-low)/3*periodictime(ipeaks)';
% adjusted_high = menvper;


% ipeaks = ipeaks(peaks > adjusted_high);

ipeaks = ipeaks(abovemed(ipeaks));


detect = zeros(size(env));


% energy considerations
% windur = 0.02;
% step = 1;
% energy = hamming_energy(xx , fs ,  windur , step);


for i=1:length(ipeaks)
    istart=ipeaks(i);
    iend = ipeaks(i);
    while(istart>1 && energy(istart)>=lowbound(istart))
        istart=istart-1;
    end
    while(iend<length(env) && energy(iend)>=lowbound(iend))
        iend = iend+1;
    end
    
    
    if(lowbound(istart) > lowbound(iend))
        tightbound = lowbound(istart);
        while(energy(iend) < tightbound)
            iend = iend-1;
        end
    else
        tightbound = lowbound(iend);
        while(energy(istart) < tightbound)
            istart = istart+1;
        end
    end
        
    if(istart~=iend)
        detect(istart:iend) = 1;
    end
end



%using area of lobe
% for i=1:length(ipeaks)
% %     disp(['--peak: ' , num2str(ipeaks(i))]);
%     istart=ipeaks(i)-1;
%     iend=ipeaks(i)+1;
%     
%     larea = env(ipeaks(i)); %lobe area 
%     flag = 1;
%     while(flag)
%         flag=0;
% %         disp(['start - ' , num2str(istart ),' : ' ,num2str(env(istart)/larea)]);
% %         disp(['end - ' , num2str(iend ),' : ' ,num2str(env(iend)/larea)]);
%         slope = env(iend+1)-env(iend);
%         if(env(iend)/larea > 10^-3 || slope < -0.01)
%             iend = iend+1;
%             larea = larea + env(iend);
%             flag=1;
%         end
%         slope = env(istart-1)-env(istart);
%         if(env(istart)/larea > 10^-3 || slope < -0.01)
%             istart = istart-1;
%             larea = larea + env(istart);
%             flag=1;
%         end
%     end
% %     figure(5)
% %     plot((istart:iend)/fs , env(istart:iend));
% %     pause(0.1)
%     
%     if(iend-istart>3)
%         detect(istart:iend)=1;
%     end
% end



%using given lower bound
%changed on 14/11/17
% for i=1:length(ipeaks)
%     istart=ipeaks(i);
%     iend = ipeaks(i);
%     while(istart>1 && env(istart)>=low(istart))
%         istart=istart-1;
%     end
%     while(iend<length(env) && env(iend)>=low(iend))
%         iend = iend+1;
%     end
%     
%     if(istart~=iend)
%         detect(istart:iend) = 1;
%     end
% end
