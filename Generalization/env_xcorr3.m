function [ periods , peaklocs , correlations ] = env_xcorr3( x , fs , env , isplot)
% Estimate periodicity (e.g. trill rate) according to normalized cross-
% correlation coefficient in the signal envelope . In case env in not a
% paramerer, it is calculated in a naive manner using the hilbert transform
%
% inputs:
% x: input signal
% fs: sample frequency
% env: signal envelope, e.g. smooth TEO or short time Hamming energy
% isplot: boolean. toggles debug plotting
%
% outputs:
% periods: array of detected periods
% peaklocs: array of respective peak locations
% correlations: array of correlation coefficient values of respective peak
% locations

maxrate = 75; %hz
minrate = 6; %hz
search_area = 0.6; %percents


if(nargin < 4)
    env = abs(hilbert(x));    
end

maxper = floor(fs / minrate);
minper = floor(fs / maxrate);


margin_env = [env(maxper:-1:1) ; env(end-maxper+1:end)];
margin_env = medfilt1(env , floor(1*fs));

%% Start from peak close enough to the "bulk of trill", calculated by using a lowpass filter with 0.5 sec time support
filtenv = filter(ones(floor(floor(fs)/2) , 1)/(floor(fs)/2) , 1 , env);
filtdelay = floor(fs/4);[~ , filtmaxt] = findpeaks(filtenv , 'SortStr', 'descend');
filtmaxt = filtmaxt-filtdelay;

[~ , t] = findpeaks(env , 'SortStr','descend');
t = t(t > maxper & t < length(env) - 1.5*maxper);
for i=1:length(t)
    if( abs( t(i) - filtmaxt(1) ) < 0.3*fs) , break; end
end

peakloc = t(i);

if(isplot)
figure(11)
plot(filtenv)
end


%% go forwards
freqparams.maxper = maxper;
freqparams.minper = minper;
freqparams.search_area = [];
[periods , peaklocs , correlations] = ...
    subfunc_env_xcorr(env , freqparams , peakloc , margin_env , isplot);

%% go backwards
backind = peaklocs(1);
[periods_back , peaklocs_back , correlations_back] = ...
    subfunc_env_xcorr(env(end:-1:1), freqparams , length(env)-backind , margin_env(end:-1:1), isplot);

peaklocs_back = length(env) - peaklocs_back; 
peaklocs = [peaklocs_back(end:-1:2) ; peaklocs];

%% 2nd iteration
freqparams.search_area = search_area;
freqparams.median_period = prctile(peaklocs(2:end)-peaklocs(1:end-1) , 70);
[periods , peaklocs , correlations] = ...
    subfunc_env_xcorr(env , freqparams, peakloc, margin_env , isplot);
backind = peaklocs(1);
[periods_back , peaklocs_back , correlations_back] = ...
    subfunc_env_xcorr(env(end:-1:1), freqparams , length(env)-backind , margin_env(end:-1:1) , isplot);

%% pack everything up
peaklocs_back = length(env) - peaklocs_back; 
periods = [periods_back(end:-1:1) ; periods];
peaklocs = [peaklocs_back(end:-1:2) ; peaklocs];
correlations = [correlations_back(end:-1:1) ; correlations];
periods = periods / fs;
peaklocs = peaklocs / fs;


function [periods , peaklocs , normcorrs] = subfunc_env_xcorr(env , freqparams , initpeakloc , margin_env, isplot)
% subfunction for normalized cross correlation (STECC) coefficients
% calculation and sequential processing
% inputs:
% env: smooth envelope of reference signal
% freqparams: structure of frequency parameters {maxper, minper,
% search_area, median_period}
% initpeakloc: initial peak location. Should be found as the maximal peak
% in the maximal energy distributed area of the signal's time
% representation.
% margin_env: envelope of the signal margins. Usually by concatenating the 
% maximal period from both signal ends.
% isplot: boolean. toggles debug plotting
%
% outputs:
% periods: estimated periods of input signal syllables
% peaklocs: peak locations of input signal syllables
% normcorrs: normalized correlation coefficients (STECC) of syllable
% centerd at peaks


margin_threshold = 1.2; % peak ratio thresold for subfuc stopping criteria

if(isempty(freqparams.search_area))
    mmaxper = freqparams.maxper;
    mminper = freqparams.minper;
else
    search_area = freqparams.search_area;
    median_period = freqparams.median_period;
    mminper = floor( (1-search_area)*median_period );
    mmaxper = floor( (1+search_area)*median_period );
end

periods = [];
peaklocs = [initpeakloc];
normcorrs = [];
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
        fprintf('env_xcorr3.subfunc_env_xcorr: Finished going over all peaks\n');
    end
    
    if(isempty(t))
        break; end
    % for plotting correlated envelopes only
    per = tpeak+mminper-1;%per = t(1)+mminper-1;
    halfper = floor(per/2);
    interval1 = peakloc + ( (-halfper):halfper );
    interval2 = interval1(end) + ( 1:length(interval1) ) ;
%     figure(10000)
%     plot(env(interval1)/norm(env(interval1)))
%     hold on
%     plot(env(interval2)/norm(env(interval2))); hold off
    if(isempty(c)) 
        break; end
    c = c(1);
    if(c < 0.7)
        break ; end

    [max_env , max_env_s] = max(env(interval2));
    if(env(peaklocs(end) + per) < margin_env(peaklocs(end) + per)*margin_threshold)
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
    if(isplot)
        figure(120)
        plot(env)
        hold on
        plot(peaklocs , env(peaklocs) , '*')
        plot(mmaxper*[1 , 1] , [0 , max(env)] , 'k');
        plot(margin_env*margin_threshold , 'm');
        plot((length(env) - mmaxper)*[1, 1] , [0 , max(env)] , 'k');
        hold off
        drawnow
    end
end
