function [ periods , peaklocs , correlations ] = env_xcorr( x , fs , time , env )
% Estimate periodicity (e.g. trill rate) according to cross-correlation in
% the signal envelope. In case env in not a parameter, it is calculated.

maxrate = 75; %hz
minrate = 7; %hz
search_area = 0.4; %percents
energy_detect = 0.15; % percents

if(nargin < 4)
    env = abs(hilbert(x));    
end

maxper = floor(fs / minrate);
minper = floor(fs / maxrate);

[~ , t] = findpeaks(env , 'SortStr','descend');
peakloc = t(1);



%% acquisition

figure(10)
plot(time , env);
corrvals = zeros(maxper - minper + 1, 1);

for per=minper:maxper
    halfper = floor(per/2);
    interval1 = peakloc + ( (-halfper):halfper );
    interval2 = interval1(end) + ( 1:length(interval1) ) ;
    corrvals(per - minper + 1) = env(interval1)' * env(interval2) / ( norm(env(interval1)) * norm(env(interval2)));
end

figure(11)

findpeaks(corrvals)
[c , t] = findpeaks(corrvals , 'SortStr','descend');

initper = t(1) + minper;

halfper = floor(initper/2);
interval1 = peakloc + ( (-halfper):halfper );
interval2 = interval1(end) + ( 1:length(interval1) ) ;
c2 = env(interval1)' * env(interval2);


[periods , peaklocs , correlations] = ...
    subfunc_env_xcorr(env , search_area , initper , peakloc , c(1) , c2);
%% go forwards
% periods = initper;
% peaklocs = [peakloc ; peakloc + initper];
% correlations = c(1);
% while(peaklocs(end) < length(env) )
%     peakloc = peaklocs(end);
%     mminper = floor( (1-search_area)*periods(end) );
%     mmaxper = floor( (1+search_area)*periods(end) );
%     corrvals = zeros(mmaxper - mminper + 1, 1);
%     for per=mminper:mmaxper
%         halfper = floor(per/2);
%         interval1 = peakloc + ( (-halfper):halfper );
%         interval2 = interval1(end) + ( 1:length(interval1) ) ;
%         corrvals(per - mminper + 1) = env(interval1)' * env(interval2) / ( norm(env(interval1)) * norm(env(interval2)));
%     end
%     findpeaks(corrvals , 'SortStr','descend')
%     [c , t] = findpeaks(corrvals , 'SortStr','descend');
%     if(isempty(c)) , break; end
%     c = c(1);
%     if(c < 0.7 || c < 0.7*mean(correlations) ) , break ; end
%     peakval = env(peaklocs(end) + t(1)+mminper);
%     meanpeak = mean(env(peaklocs));
% %     if(peakval < energy_detect * meanpeak), break; end
%     periods = [periods; t(1)+mminper];
%     peaklocs = [peaklocs ; peaklocs(end) + periods(end)];
%     correlations = [correlations; c];
% end

backind = peaklocs(2);
[periods_back , peaklocs_back , correlations_back] = ...
    subfunc_env_xcorr(env(end:-1:1), search_area , initper , length(env)-backind , c(1));

peaklocs_back = 2*peaklocs_back(1) - peaklocs_back; 
% periods_back = periods(1);
% peaklocs_back = peaklocs(1);
% correlations_back = correlations(1);
% %% go backwards
% periods_back = periods(1);
% peaklocs_back = peaklocs(1);
% correlations_back = correlations(1);
% while(peaklocs_back(end) > 1)
%     peakloc = peaklocs_back(end);
%     mminper = floor( (1-search_area)*periods_back(end) );
%     mmaxper = floor( (1+search_area)*periods_back(end) );
%     corrvals = zeros(mmaxper - mminper + 1, 1);
%     for per=mminper:mmaxper
%         halfper = floor(per/2);
%         interval1 = peakloc + ( (-halfper):halfper );
%         interval2 = interval1(1) - ( 1:length(interval1) ) ;
% %         norm(env(interval1)) / norm(env(interval2))
%         corrvals(per - mminper + 1) = env(interval1)' * env(interval2) / ( norm(env(interval1)) * norm(env(interval2)));
%     end
%     findpeaks(corrvals , 'SortStr','descend')
%     [c , t] = findpeaks(corrvals , 'SortStr','descend');
%     if(isempty(c)) , break; end
%     c = c(1);
%     if(c < 0.7 || c < 0.7*mean(correlations_back) ) , break ; end
%     peakval = env(peaklocs_back(end) - ( t(1)+mminper ) );
%     meanpeak = mean(env(peaklocs_back));
% %     if(peakval < energy_detect * meanpeak) , break; end
% %     env(peaklocs_back(end)) / env(peaklocs_back(end)+-(t(1)+mminper))
%     periods_back = [periods_back; t(1)+mminper];
%     peaklocs_back = [peaklocs_back ; peaklocs_back(end) - periods_back(end)];
%     correlations_back = [correlations_back; c];
%     
% end
periods = [periods_back(end:-1:2) ; periods];
peaklocs = [peaklocs_back(end:-1:2) ; peaklocs];
correlations = [correlations_back(end:-1:2) ; correlations];

periods = periods / fs;
peaklocs = peaklocs / fs;


function [periods , peaklocs , normcorrs] = subfunc_env_xcorr(env , search_area , initper , initpeakloc , initnormcorrs , initcorr2)
% subfunction for self cross correlation
periods = initper;
peaklocs = [initpeakloc ; initpeakloc + initper];
normcorrs = initnormcorrs;
correlations2 = initcorr2;
while(peaklocs(end) < length(env) )
    peakloc = peaklocs(end);
    mminper = floor( (1-search_area)*periods(end) );
    mmaxper = floor( (1+search_area)*periods(end) );
    corrvals = zeros(mmaxper - mminper + 1, 1);
    for per=mminper:mmaxper
        halfper = floor(per/2);
        interval1 = peakloc + ( (-halfper):halfper );
        interval2 = interval1(end) + ( 1:length(interval1) ) ;
        corrvals(per - mminper + 1) = env(interval1)' * env(interval2) / ( norm(env(interval1)) * norm(env(interval2)));
    end
%     findpeaks(corrvals , 'SortStr','descend')
    [c , t] = findpeaks(corrvals , 'SortStr','descend');
    
    % for plotting correlated envelopes only
    per = t(1)+mminper-1;
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
    if(c < 0.7 || c < 0.7*mean(normcorrs) )
        break ; end
    peakval = env(peaklocs(end) + t(1)+mminper);
    meanpeak = mean(env(peaklocs));
%     if(peakval < energy_detect * meanpeak), break; end
    periods = [periods; t(1)+mminper];
    peaklocs = [peaklocs ; peaklocs(end) + periods(end)];
    normcorrs = [normcorrs; c];
    correlations2 = [correlations2; c2];
    figure(12)
    plot(env)
    hold on
    plot(peaklocs , env(peaklocs) , '*')
    hold off
    drawnow
end
