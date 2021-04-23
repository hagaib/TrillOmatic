function [ detect , segs , env] = longtrill_syllable_detectionAOI( x , fs , yin , aoipeaks , f0band , isplot)
%LONGTRILL_SYLLABLE_DETECTIONAOI Summary of this function goes here
%   Detailed explanation goes here

%to account for delay bug in double filtered xteo in area_of_interest4
% delay = 0.02;
% aoipeaks = aoipeaks - delay;
time = (0:length(x)-1)/fs;

periods = aoipeaks(2:end) - aoipeaks(1:end-1);

%% TEO estimation
wpass = f0band /fs *2;
wstop = [wpass(1)*0.98 , wpass(2)*1.02];
ftype = 'bandpass';
if(wpass(2) >=1 || wstop(2) >=1 ) %high pass
    wpass = wpass(1); wstop = wstop(1);
    ftype = 'high';
end
Rp = 0.05; Rs = 40;
n = cheb2ord(wpass , wstop , Rp , Rs);
[z,p,g] = cheby2(n,Rs , wstop , ftype);
sos = zp2sos(z,p,g);
xx0 = sosfilt(sos , x);

energy = xx0.^2;
winlen = median(periods);
s_winlen = floor(winlen/2*fs); s_winlen = s_winlen + 1-mod(s_winlen , 2);
s_winlenShort = floor(s_winlen/4) + 1 - mod(floor(s_winlen/4) , 2);
hammingEnergy = filter(hamming(s_winlen)/sum(hamming(s_winlen)) , 1 , energy);
hammingEnergyShort = filter(hamming(s_winlenShort)/sum(hamming(s_winlenShort)) , 1 , energy);
hammingEnergySteadystate = hammingEnergy(ceil(s_winlen/2):end);
hammingEnergy = [hammingEnergy(ceil(s_winlen/2):end) ; zeros(floor(s_winlen/2) , 1)];
hammingEnergyShort = [hammingEnergyShort(ceil(s_winlenShort/2):end) ; zeros(floor(s_winlenShort/2) , 1)];

if(isplot)
figure(1010)
plot(energy , 'b'); hold on; plot(hammingEnergy , 'k');
plot(hammingEnergyShort , 'g');
plot([1,length(energy)],1.5*median(hammingEnergySteadystate)*[1,1],'r'); 
plot([1,length(energy)],1.5*median(hammingEnergyShort)*[1,1],'m'); hold off
end
% [~ , threshold_teo , teotime] = trilltime_TEO(xx0 , fs);

teo = teager(xx0);
teo = [teo(1) ; teo ; teo(end)];
L = floor(fs*0.02);
if(mod(L,2)==0) , L=L+1; end
g = gausswin(L);
g = g/sum(g);
teo = filter(g,1,teo);
teo = filter(g,1,teo(end:-1:1)); teo = teo(end:-1:1);

newNyquist = fs/s_winlen*30;
subsampleFactor = floor(fs/(newNyquist*2));
hammingEnergySubsample = hammingEnergy(1:subsampleFactor:end);

[periods , peaklocs] = env_xcorr3([] , fs/subsampleFactor, time(1:subsampleFactor:end), hammingEnergySubsample , 5, isplot);

zcrwindur = 1/f0band(1)*5; % 5 cycles of the longest possible period
winsamples = floor(-zcrwindur*fs/2):floor(zcrwindur*fs/2);
winsamples_yin = floor(-zcrwindur*yin.steprate/2):floor(zcrwindur*yin.steprate/2);
peak_zcr = zeros(size(peaklocs));
peak_yin = zeros(size(peaklocs));
peak_eng = zeros(size(peaklocs));
for i=1:length(peak_zcr)
    samples = floor(fs*peaklocs(i)) + winsamples;
    peak_zcr(i) = zcr(x(samples),1/fs,1);
    peak_yin(i) = median(yin.f0(find(yin.time>peaklocs(i) , 1)+winsamples_yin),  'omitnan');
    peak_eng(i) = max(hammingEnergy(samples));
end

% %% find ZCR outliers using quartiles
q25 = prctile(peak_zcr,25); q75 = prctile(peak_zcr,75);
% iqr_zcr = q75 - q25;
% cutoff = iqr_zcr * 3;
% lower = q25 - cutoff; upper = q75 + cutoff;
% 
% outliers = peak_zcr < lower | peak_zcr > upper;


% % maybe with normal distribution approximation?
% median_yin = median(peak_yin(floor(end/4):floor(end*3/4)));
% stdev = std(peak_yin);
% cutoff_stds = 3;
% lower = median_yin - cutoff_stds*stdev; upper = median_yin + cutoff_stds*stdev;


%% last chance: outliers with zcr at margins


margin_t = 0.1; %seconds
margin_s = floor(margin_t * fs);
margin_zcr_l = zcr(x(1:1+margin_s),1/fs,1);
margin_zcr_r = zcr(x(end-margin_s+1:end),1/fs,1);
% margin_zcr_mean = mean([margin_zcr_l , margin_zcr_r]);
start_i=1; end_i=length(peaklocs);

if(margin_zcr_l < q25 || margin_zcr_l > q75)
    for start_i=1:length(peaklocs)
        if( ~( ((peak_zcr(start_i) < margin_zcr_l && margin_zcr_l < q25)||...
                (peak_zcr(start_i) > margin_zcr_l && margin_zcr_l > q75)  ) &&...
                peak_eng(start_i) < median(peak_eng)/10 ))
            break;
        end
%         if( ((peak_zcr(start_i) > margin_zcr_l || margin_zcr_l > q75) && ...
%             (peak_zcr(start_i) < margin_zcr_l || margin_zcr_l < q25)) || ...
%                 peak_eng(start_i)>median(peak_eng)/10) , break; end   
    end
end
if(margin_zcr_r < q25 || margin_zcr_r > q75)
    for end_i=length(peaklocs):-1:start_i
        if( ~( ((peak_zcr(end_i) < margin_zcr_r && margin_zcr_r < q25)||...
                (peak_zcr(end_i) > margin_zcr_r && margin_zcr_r > q75)  ) &&...
                peak_eng(end_i) < median(peak_eng)/10 ))
            break;
        end
%         if(peak_zcr(end_i) < margin_zcr_r || peak_eng(end_i)>median(peak_eng)/10) , break; end    
    end
end
outliers = 1:length(peaklocs) < start_i | 1:length(peaklocs) > end_i;


if( sum(outliers) )
%     if(isplot)
        fprintf('Peak at %0.3f was recognized as outlier and removed\n' , peaklocs(outliers));
%     end
    peaklocs = peaklocs(~outliers);
end


[detect , segs , plotparams] = longtrill_segmentation(x , fs , hammingEnergy, peaklocs , 'adaptive' , hammingEnergyShort);

% %% Detecting pulse interval
% detect = zeros(size(x));
% windur = 0.05;
% % low = prcfilt(teotime, floor(windur*fs) , 5 , 'omitnan')+threshold_teo;
% low = 2*median(hammingEnergyShort);
% 
% for i=1:length(peaklocs)
%     istart = floor(peaklocs(i)*fs);
%     iend = floor(peaklocs(i)*fs);
%     while(istart>0 && hammingEnergyShort(istart)>low) , istart = istart-1; end
%     while(iend<length(x) && hammingEnergyShort(iend)>low) , iend = iend+1; end
%         
%     detect(istart:iend) = 1;
% end
% 
% detect = medfilt1(detect , 3);
% segs = logical2segments(detect , fs);

if(isplot)
figure(301)
plot(time , xx0*10^-1 , 'b'); hold on;
plot(time , hammingEnergyShort , 'k-');
plot(time , hammingEnergy , 'k.');
plot(time , plotparams.low , 'g'); 
plot(time , detect*max(teo) , 'm') , 
plot(peaklocs , teo(floor(peaklocs*fs)) , '*'); hold off;
end

env = hammingEnergyShort;

end

