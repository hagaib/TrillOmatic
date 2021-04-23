function [ detect , segs , env ] = longtrill_syllable_detection2(xx , fs , yin , harmonics , harmenergy , axes , f0band)
%%
%this version focuses on energy envelope methods
%%

%% 2 Options here:
% 1) Let trilltime estimation be strict , locating only the part where the
% longtrill is present within the recording. Drawback: Usually fails to
% find weak pulses. Ocassionally misses more.
% 2) Let trilltime estimation be loose. Afterwards detect every pulse in
% recording , regardless of connection to longtrill. After initial
% detection, a rejection process based on yin periodicity and gap variance
% eliminates erroneous detections. Drawback: Variance analysis can be very
% often misleading, especially with other bird sounds present in recording

if(nargin < 7)
    f0band = [2000 , 3500];
end

time = ((0:length(xx)-1)/fs)';
min_segment_time = 0.015;
[periodictimeyin , nonperiodictimeyin , persegs , nonpersegs] = ...
    pertime(yin.dips , yin.time , min_segment_time , yin.f0 , f0band );

% [minf0 , maxf0 ,minf1 , maxf1 , minf2 , maxf2] = harmonic_boundaries(yin , harmonics , periodictimeyin , f0band);
% fbounds = [minf0 , maxf0; minf1 , maxf1 ; minf2, maxf2];

xx0 = filter_signal(xx , fs , f0band(1) , f0band(2));

energy = hamming_energy(xx0, fs ,  0.01 , 1);

bulktime = longtrill_bulk_of_mass(persegs , yin , 0.5 , 1)
% [istart , iend] = deal(bulktime(1) , bulktime(2));

trilltime = trilltime_estimate(xx0 , fs , yin , 0.01 , 200)
[tstart , tend] = deal(min(trilltime(1) , bulktime(1)) , max(trilltime(2) , bulktime(2)));
%% For plots made for article only
% tstart = tstart - 0.05;
% tend = tend + 0.1;


% outoftrill = nonpersegs<yin.time(istart) | nonpersegs>yin.time(iend);
% nonpertime = segments2logical(nonpersegs(:,outoftrill(1,:)) , time , fs);
% energyn = mean(energy(nonpertime));

[env , detect , segs , PRI] = energy_envelope2(energy , fs , time , axes , [tstart , tend] , xx0);
% [tstart , tend] = longtrill_fine_time2(fs, yin , segs , energy , energyn , PRI);
% [tstart , tend] = logtrill_fine_time(fs , yin , segs , energy , energyn, env);
% [env , detect , segs] = energy_envelope2(energy , fs , time , axes , [tstart , tend] , xx0);



%% YIN Periodicity Rejects
dipthresh = 0.3;
segdips = zeros(size(segs,2),1);
yin.dips(isnan(yin.dips)) = 1;
for i=1:size(segs , 2)
    seg = segs(:,i);
    segdips(i) = mean(yin.dips(yin.time>=seg(1) & yin.time <=seg(2)));
end
segs = segs(:,segdips<dipthresh);


%% Length Rejects
lenthresh = 0.008;
seglen = segs(2,:) - segs(1,:);
segs = segs(: , seglen>lenthresh);


% %% Variance Rejects
% diptime = mean(segs);
% gaptime = diptime(2:end) - diptime(1:end-1);
% [mg,ig] = max(gaptime);
% while mg > PRI + 4*std(gaptime)
%     center = find(diptime > mean(diptime) , 1);
%     if(ig>center)
%         ig = ig+1;
%     end
%     segs(:,ig) = [];
%     diptime(ig) = [];
%     gaptime = diptime(2:end) - diptime(1:end-1);
%     [mg,ig] = max(gaptime);
% end


function [tstart , tend] = longtrill_fine_time2(fs , yin , detectsegs , energy , noiseenergy , PRI)
%same as longtrill_fine_time but considering new syllables based on yin and
%energy data and not on persegs.

energythresh = noiseenergy;

yinfs = 1/(yin.time(2) - yin.time(1));

% PRI = median(detectsegs(1,2:end) - detectsegs(1 , 1:end-1));
% mdiff = mean(PRI);
syllength = detectsegs(2,:) - detectsegs(1,:);

tstart = detectsegs(1,1);
tend = detectsegs(2 , end);

memory = 3;
pulsestart = tstart;
meanyin = mean(yin.dips(floor(yinfs*detectsegs(1,1)):floor(yinfs*detectsegs(1,1+memory))) , 'omitnan');
meaneng = mean(energy(floor(fs*detectsegs(1,1)):floor(fs*detectsegs(1,1+memory))));
pri = PRI;
while(pulsestart(1)-pri*1.1>0)
    pri = pri*1.1; % pulse spacings tend to accelerate at the begining
    newpulsetime = [pulsestart(1)-pri , pulsestart(1)];
    newyin = mean(yin.dips(floor(yinfs*(pulsestart(1)-pri)):floor(yinfs*pulsestart(1))) , 'omitnan');
    neweng = mean(energy(floor(fs*newpulsetime(1)):floor(fs*newpulsetime(2))));
    if(newyin < meanyin*1.8 || neweng > meaneng*0.3)
        pulsestart = [newpulsetime(1) , pulsestart(1)];
        detectsegs = [pulsestart' , detectsegs];
        meanyin = mean(yin.dips(floor(yinfs*detectsegs(1,1)):floor(yinfs*detectsegs(1,1+memory))) , 'omitnan');
        meaneng = mean(energy(floor(fs*detectsegs(1,1)):floor(fs*detectsegs(1,1+memory))));
    else
        break
    end
end
tstart = pulsestart(1);

total_time = (length(energy)-1)/fs;
pulseend = tend;
meanyin = mean(yin.dips(floor(yinfs*detectsegs(2,end-memory)):floor(yinfs*detectsegs(2,end))));
meaneng = mean(energy(floor(fs*detectsegs(2,end-memory)):floor(fs*detectsegs(2,end))));
pri = PRI;
while(pulseend(end)+pri<total_time)
    pri = pri*1.1; % pulse spacings tend to accelerate at the begining
    newpulsetime = [pulseend(end), pulseend(end)+pri];
    newyin = mean(yin.dips(floor(yinfs*pulseend(end)):floor(yinfs*(pulseend(end)+pri))));
    neweng = mean(energy(floor(fs*newpulsetime(1)):floor(fs*newpulsetime(2))));
    if(newyin < meanyin*1.6 || neweng > meaneng*0.3)
        pulseend = [pulseend(end) , newpulsetime(2)];
        detectsegs = [detectsegs, pulseend'];
        meanyin = mean(yin.dips(floor(yinfs*detectsegs(2,end-memory)):floor(yinfs*detectsegs(2,end))));
        meaneng = mean(energy(floor(fs*detectsegs(2,end-memory)):floor(fs*detectsegs(2,end))));
    else
        break
    end
end
tend = pulseend(end);



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
n = cheb2ord([flow , fhigh]/fs*2 , [flow*0.95 , fhigh*1.05]/fs*2 , 0.5 , 40);
[z,p,k] = cheby2(n , 40, [flow , fhigh]/fs*2 , 'bandpass');
sos = zp2sos(z,p,k);
xxfilt = sosfilt(sos , xx);

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
