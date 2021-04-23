function [ detect , segs , env  , smoothenergy] = longtrill_syllable_detection3(xx , fs , yin , harmonics , harmenergy , axes)
%%
%this version focuses on energy envelope methods
%%
time = ((0:length(xx)-1)/fs)';
min_segment_time = 0.015;
[periodictimeyin , nonperiodictimeyin , persegs , nonpersegs] = ...
    pertime(yin.dips , yin.time , min_segment_time , yin.f0 , [2000 , 3500] );

[minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2] = harmonic_boundaries(yin , harmonics , periodictimeyin , [2000 , 3500]);
fbounds = [minf0 , maxf0; minf1 , maxf1 ; minf2, maxf2];

xx0 = filter_signal(xx , fs , fbounds(1,1) , fbounds(1,2));

energy = hamming_energy(xx0, fs ,  0.01 , 1);
teo = teager(xx0);
teo = [teo(1);teo;teo(end)];

bulktime = longtrill_bulk_of_mass(persegs , yin , 1.5 , 0);
[istart , iend] = deal(bulktime(1) , bulktime(2));
 
 
outoftrill = nonpersegs<yin.time(istart) | nonpersegs>yin.time(iend);
nonpertime = segments2logical(nonpersegs(:,outoftrill(1,:)) , time , fs);
energyn = mean(energy(nonpertime));


[env , detect , segs] = energy_envelope3(teo , fs , time , [] , yin.time([istart , iend]) , xx0);
[tstart , tend] = longtrill_fine_time2(fs, yin , segs , energy , energyn);
% [tstart , tend] = logtrill_fine_time(fs , yin , segs , energy , energyn, env);
[env , detect , segs , smoothenergy] = energy_envelope3(teo , fs , time , axes , [tstart , tend] , xx0);




function [tstart , tend] = longtrill_fine_time2(fs , yin , detectsegs , energy , noiseenergy)
%same as longtrill_fine_time but considering new syllables based on yin and
%energy data and not on persegs.

energythresh = noiseenergy;

yinfs = 1/(yin.time(2) - yin.time(1));

PRI = median(detectsegs(1,2:end) - detectsegs(1 , 1:end-1));
mdiff = mean(PRI);
syllength = detectsegs(2,:) - detectsegs(1,:);

tstart = detectsegs(1,1);
tend = detectsegs(2 , end);

memory = 3;
pulsestart = tstart;
meanyin = mean(yin.dips(floor(yinfs*detectsegs(1,1)):floor(yinfs*detectsegs(1,1+memory))));
meaneng = mean(energy(floor(fs*detectsegs(1,1)):floor(fs*detectsegs(1,1+memory))));
while(pulsestart(1)-PRI>0)
    newpulsetime = [pulsestart(1)-PRI , pulsestart(1)];
    newyin = mean(yin.dips(floor(yinfs*(pulsestart(1)-PRI)):floor(yinfs*pulsestart(1))));
    neweng = mean(energy(floor(fs*newpulsetime(1)):floor(fs*newpulsetime(2))));
    if(newyin < meanyin*1.2 || neweng > meaneng*0.8)
        pulsestart = [newpulsetime(1) , pulsestart];
    else
        break
    end
end
tstart = pulsestart(1);


function [tstart , tend] = logtrill_fine_time(fs , yin , detectsegs , energy , noiseenergy , env)
energythresh = noiseenergy;

syldiff = detectsegs(1,2:end) - detectsegs(1 , 1:end-1);
mdiff = max(syldiff);
syllength = detectsegs(2,:) - detectsegs(1,:);

[periodictimeyin , ~ , persegs , ~] = pertime(yin.dips , yin.time , 0.005 , yin.f0 ,  [2000 , 3500] , 0.3);

outsegs = find(persegs(2,:) >= detectsegs(2,end));
lastseg = detectsegs(:,end);
for i=1:length(outsegs)
    if( persegs(1,outsegs(i)) - lastseg(1) > mdiff * 1.5) , break; end
    yinind = yin.time>persegs(1,outsegs(i)) & yin.time<persegs(2,outsegs(i));
    xxind = floor(persegs(:,outsegs(i))*fs);
    dips = yin.dips(yinind);
    eng = energy(xxind(1):xxind(2));
    
    if(mean(eng) > energythresh)
        lastseg = persegs(: , outsegs(i));
    end
    
end
% tend = find(yin.time > lastseg(2) , 1 , 'first');
tend = lastseg(2);


outsegs = find(persegs(1,:) <= detectsegs(1,1));
lastseg = detectsegs(:,1);
for i=length(outsegs):-1:1
    if( lastseg(1) - persegs(1,outsegs(i)) > mdiff * 1.5) , break; end
    yinind = yin.time>persegs(1,outsegs(i)) & yin.time<persegs(2,outsegs(i));
    xxind = floor(persegs(:,outsegs(i))*fs);
    dips = yin.dips(yinind);
    eng = energy(xxind(1):xxind(2));
    
    if(mean(eng) > energythresh)
        lastseg = persegs(: , outsegs(i));
    end
    
end
% tstart = find(yin.time > lastseg(1) , 1 , 'first');
tstart = lastseg(1);




function xxfilt = filter_signal(xx , fs , flow , fhigh)
%filtering the signal to desired spectrum band
wn = [flow , fhigh]/(fs/2);
wn(end) = min(wn(end) , 0.9);
[b,a] = butter(5 , wn , 'bandpass');
% [b , ~] = fir1(72 , wn , 'bandpass'); a=1;
xxfilt = filtfilt(b , a , xx);

function [minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2] = harmonic_boundaries(yin , harmonics , periodictimeyin , f0boundaries)


[minf0 , maxf0] = bw_boundaries (yin.f0 , periodictimeyin, 800, 20);
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

function [minf , maxf] = bw_boundaries (fvect , selection , maxbw , p)
%fvect = vector of frequencies
%selection = vector of 1's and 0's (indices chosen)
%p = (optional) precentile above minimum 

if (nargin < 4) , p=0; end

perf = fvect(selection);
minf = prctile(perf , p);
maxf = prctile(perf , 100-p);

[minf , maxf] = clip_to_min_bw(minf , maxf , maxbw);
