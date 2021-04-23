function [ periods , peaklocs , correlations ] = env_xcorr2( x , fs , time , env , corr2threshold, isplot)
% Estimate periodicity (e.g. trill rate) according to cross-correlation in
% the signal envelope. In case env in not a parameter, it is calculated.

maxrate = 75; %hz
minrate = 6; %hz
search_area = 0.6; %percents

if(nargin < 4)
    env = abs(hilbert(x));    
end

maxper = floor(fs / minrate);
minper = floor(fs / maxrate);


%% Start from peak close enough to the "bulk of trill", calculated by using a lowpass filter with 0.5 sec time support
filtenv = filter(ones(floor(floor(fs)/2) , 1)/(floor(fs)/2) , 1 , env);
[~ , filtmaxt] = findpeaks(filtenv , 'SortStr', 'descend');

[~ , t] = findpeaks(env , 'SortStr','descend');
t = t(t > maxper & t < length(env) - 1.5*maxper);
for i=1:length(t)
    if( abs( t(i) - filtmaxt(1) ) < 0.3*fs) , break; end
end

peakloc = t(i);

if(isplot)
figure(11)
plot(time , filtenv)
end

%% acquisition
if(isplot)
figure(10)
plot(time , env);
end
corrvals = zeros(maxper - minper + 1, 1);

for per=minper:maxper
    halfper = floor(per/2);
    interval1 = peakloc + ( (-halfper):halfper );
    interval2 = interval1(end) + ( 1:length(interval1) ) ;
    if(interval2(end) > length(env))
        corrvals = corrvals(1:per - minper);
        break;
    end
    corrvals(per - minper + 1) = env(interval1)' * env(interval2) / ( norm(env(interval1)) * norm(env(interval2)));
end

% figure(11)
% findpeaks(corrvals)
try
[c , t] = findpeaks(corrvals , 'SortStr','descend');
catch
    fprintf('yo mama\n');
end

if(length(t)<1)
    [ periods , peaklocs , correlations ] = deal([],[],[]);
    return;
end
initper = t(1) + minper;

halfper = floor(initper/2);
interval1 = peakloc + ( (-halfper):halfper );
interval2 = interval1(end) + ( 1:length(interval1) ) ;
c2 = env(interval1)' * env(interval2);

%% go forwards
freqparams.maxper = maxper;
freqparams.minper = minper;
freqparams.search_area = [];
[periods , peaklocs , correlations] = ...
    subfunc_env_xcorr(env , freqparams , initper , peakloc , c(1) , c2 , corr2threshold , isplot);

%% go backwards
backind = peaklocs(2);
[periods_back , peaklocs_back , correlations_back] = ...
    subfunc_env_xcorr(env(end:-1:1), freqparams , initper , length(env)-backind , c(1) , c2, corr2threshold, isplot);

peaklocs_back = length(env) - peaklocs_back; 
peaklocs = [peaklocs_back(end:-1:3) ; peaklocs];

%% 2nd iteration
freqparams.search_area = search_area;
freqparams.median_period = prctile(peaklocs(2:end)-peaklocs(1:end-1) , 70);
[periods , peaklocs , correlations] = ...
    subfunc_env_xcorr(env , freqparams , initper , peakloc , c(1) , c2, corr2threshold , isplot);
backind = peaklocs(2);
[periods_back , peaklocs_back , correlations_back] = ...
    subfunc_env_xcorr(env(end:-1:1), freqparams , initper , length(env)-backind , c(1) , c2, corr2threshold , isplot);

%% pack everything up
peaklocs_back = length(env) - peaklocs_back; 
periods = [periods_back(end:-1:2) ; periods];
peaklocs = [peaklocs_back(end:-1:3) ; peaklocs];
correlations = [correlations_back(end:-1:2) ; correlations];
periods = periods / fs;
peaklocs = peaklocs / fs;


function [periods , peaklocs , normcorrs] = subfunc_env_xcorr(env , freqparams , initper , initpeakloc , initnormcorrs , initcorr2 , corr2threshold, isplot)
% subfunction for self cross correlation

if(nargin<7)
    corr2threshold = 10^-1;
end
if(isempty(freqparams.search_area))
    mmaxper = freqparams.maxper;
    mminper = freqparams.minper;
else
    search_area = freqparams.search_area;
    median_period = freqparams.median_period;
    mminper = floor( (1-search_area)*median_period );
    mmaxper = floor( (1+search_area)*median_period );
end

periods = initper;
peaklocs = [initpeakloc ; initpeakloc + initper];
normcorrs = initnormcorrs;
correlations2 = initcorr2;
while(peaklocs(end) + ceil(1.5*mmaxper) < length(env) )
    peakloc = peaklocs(end);

    corrvals = zeros(mmaxper - mminper + 1, 1);
    for per=mminper:mmaxper
        halfper = floor(per/2);
        interval1 = peakloc + ( (-halfper):halfper );
        interval2 = interval1(end) + ( 1:length(interval1) ) ;
        if(interval2(end) > length(env) || interval1(1) < 1)
            corrvals = corrvals(1:per-mminper); break; end
        corrvals(per - mminper + 1) = env(interval1)' * env(interval2) / ( norm(env(interval1)) * norm(env(interval2)));
    end
    
    if(isplot)
        figure(10000)
        findpeaks(corrvals , 'SortStr','descend')
    end
    try
    [c , t] = findpeaks(corrvals , 'SortStr','descend');
    [~ , tdips] = findpeaks(-corrvals);
    tpeak = t(1);
    i=1;
    while(t(i) < tdips(1)) , i = i+1; end
    tpeak = t(min(i ,length(t)));
    catch
        fprintf('Hey what''s going on here??\n');
    end
    
    if(isempty(t))
        break; end
    % for plotting correlated envelopes only
    per = tpeak+mminper-1;%per = t(1)+mminper-1;
    halfper = floor(per/2);
    interval1 = peakloc + ( (-halfper):halfper );
    interval2 = interval1(end) + ( 1:length(interval1) ) ;
    c2 = env(interval1)' * env(interval2);
%     figure(10000)
%     plot(env(interval1)/norm(env(interval1)))
%     hold on
%     plot(env(interval2)/norm(env(interval2))); hold off
    if(isempty(c)) 
        break; end
    c = c(1);
    if(c < 0.7)% || c < 0.7*mean(normcorrs) )
        break ; end
    if(c2 < mean(correlations2) * corr2threshold)
        break; end


%     peakval = env(peaklocs(end) + t(1)+mminper);
%     meanpeak = mean(env(peaklocs));
%     if(peakval < energy_detect * meanpeak), break; end
    peakloc = peaklocs(end) + per;
    maxwidth = 0.03;
    max_search_area = peakloc + floor(-(per)*maxwidth:(per)*maxwidth);
    [~ , pp] = max(env(max_search_area));
    peakloc = peakloc + floor(-(per)*maxwidth) + pp(1);
    periods = [periods; per];
    peaklocs = [peaklocs ; peakloc];
    normcorrs = [normcorrs; c];
    correlations2 = [correlations2; c2];
    if(isplot)
        figure(120)
        plot(env)
        hold on
        plot(peaklocs , env(peaklocs) , '*')
        hold off
        drawnow
    end
end
